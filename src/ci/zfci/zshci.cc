#include <src/ci/zfci/zshci.h>
#include <src/util/math/comb.h>
#include <src/util/exception.h>
#include <src/ci/zfci/reljop.h>
#include <stdio.h>
#include <src/util/prim_op.h>
#include <src/mat1e/rel/reloverlap.h>

using namespace bagel;
using namespace std;
using namespace btas;

ZSHCI::ZSHCI(shared_ptr<const PTree> p, shared_ptr<const Geometry> g, shared_ptr<const Reference> b,
               const int ncore, const int nocc, shared_ptr<const ZCoeff_Block> coeff_zcas, const bool store_c, const bool store_g)
 : ZHarrison(p, g, b, ncore, nocc, coeff_zcas, store_c, store_g) {
  auto rr = dynamic_pointer_cast<const RelReference>(ref_);
  if (!rr) throw runtime_error("ZSHCI requires a relativistic reference boject");

  gaunt_ = idata_->get<bool>("gaunt", rr->gaunt());
  breit_ = idata_->get<bool>("breit", rr->breit());

  cout << "    * Relativistic FCI" << endl;
  cout << "    * " << nele_ << " active electrons in " << norb_ << " orbitals."  << endl;
  cout << "    * gaunt    : " << (gaunt_ ? "true" : "false") << endl;
  cout << "    * breit    : " << (breit_ ? "true" : "false") << endl;

  if (!geom_->dfs())
    geom_ = geom_->relativistic(gaunt_);

  shared_ptr<const ZCoeff_Block> coeff = coeff_zcas;
  if (!coeff)
    coeff = init_coeff();
  update(coeff);

  if (idata_->get<bool>("only_ints", false)) {
    dump_integrals();
    throw Termination("Relativistic MO are dumped onto integrals");
  }

  shci_data_ = idata_->get_child_optional("shci");
  if (!shci_data_) throw runtime_error("ZSHCI requires a subblock to initialize Dice parameters");
  if (idata_->get<string>("title") == "zcasscf") do_rdm_ = true;
  if (idata_->get<string>("title") == "zshci"  ) do_rdm_ = false;
  tight_ = false;
  first_iter_ = false;

  // A few lines to deal with Dice parameters
  // Currently I would suppose input is provided and the only option is whether to generate guess determinant for Dice input
  //generate_guess_ = idata_->get<bool>("generate_guess", false);
  rdm_file_  = shci_data_->get<string>("rdm_file", "Dice");
  dice_dir_  = shci_data_->get<string>("Dice_exe", "ZDice2");
  mpiprefix_ = shci_data_->get<string>("mpiprefix", "");
}

void ZSHCI::read_dice_rdm() {
  
}

void ZSHCI::disk_energy() {
  FILE* f = fopen("./shci.e", "rb");
  double E[10];
  fread(E, sizeof(double), nstate_, f);
  fclose(f);
  //ifstream fs("./shci.e");
  //string line;
  //while(getline(fs, line)) {
  //  stringstream ss(line);
  //  double energy;
  //  ss >> energy;
  //  energy_.push_back(energy);
  //} 
  for (int i = 0; i < nstate_; i++) {
    energy_[i] = E[i];
    cout << "state " << i << " : " << energy_[i] << endl;
  }
}

void ZSHCI::dump_integrals() const {
#ifndef DISABLE_SERIALIZATION
  auto rr = dynamic_pointer_cast<const RelReference>(ref_);
  assert(rr);
  OArchive ar("relref");
  auto rout = make_shared<RelReference>(geom_, jop_->coeff()->striped_format(), 0.0, rr->nneg(), rr->nclosed(), rr->nact(), rr->nvirt(), rr->gaunt(), rr->breit(), rr->kramers());
  ar << rout;
  dump_ints();
#else
  throw runtime_error("You must compile with serialization in order to dump MO integrals into a file.");
#endif
}

void ZSHCI::dump_integrals_and_exit() const {
  return;
}

void ZSHCI::update(shared_ptr<const ZCoeff_Block> coeff) {
  Timer timer;
  jop_ = make_shared<RelJop>(geom_, ncore_*2, (ncore_+norb_)*2, coeff, gaunt_, breit_, store_half_ints_, store_gaunt_half_ints_);
  cout << "    * Integral transformation done. Elapsed time: " << setprecision(2) << timer.tick() << endl << endl;
}

void ZSHCI::clean() {
  // cleaning up Dice temp files
  system("rm Dice_*rdm*");
  system("rm FCIDUMP");

  system("rm *.bkp");
  return;
}

shared_ptr<const ZCoeff_Block> ZSHCI::init_coeff() {
  auto rr = dynamic_pointer_cast<const RelReference>(ref_);
  assert(rr);
  auto scoeff = make_shared<const ZCoeff_Striped>(*rr->relcoeff_full(), ncore_, norb_, rr->relcoeff_full()->mdim()/4-ncore_-norb_, rr->relcoeff_full()->mdim()/2);
  
  const shared_ptr<const PTree> iactive = idata_->get_child_optional("active");
  if (iactive) {
    set<int> active_indices;
    for (auto& i : *iactive)
      active_indices.insert(lexical_cast<int>(i->data()) - 1);
    scoeff = scoeff->set_active(active_indices, geom_->nele()-charge_);
  }
  return scoeff->block_format();
}

void ZSHCI::compute() {
  Timer pdebug(2);

  /*if (!restarted_) {
    const static Comb combination;
    const size_t max_states = combination(2*norb_, nele_);
    if (nstate_ > max_states) {
      const string space = "(" + to_string(nele_) + "," + to_string(norb_) + ")";
      throw runtime_error("Wrong states specified - a " + space + " active space can only produce " + to_string(max_states) + " eigenstates.");
    }
    // find determinants that have small diagonal energies
    int offset = 0;
    for (int ispin = 0; ispin != states_.size(); ++ispin) {
      int nstate = 0;
      for (int i = ispin; i != states_.size(); ++i)
        nstate += states_[i];

      if (nstate == 0)
        continue;

      if ((geom_->nele()+ispin-charge_) % 2 == 1) {
        if (states_[ispin] == 0) {
          continue;
        } else {
          if ((geom_->nele()-charge_) % 2 == 0) throw runtime_error("Wrong states specified - only integer spins are allowed for even electron counts.");
          else throw runtime_error("Wrong states specified - only half-integer spins are allowed for odd electron counts.");
        }
      }

      const int nelea = (geom_->nele()+ispin-charge_)/2 - ncore_;
      const int neleb = (geom_->nele()-ispin-charge_)/2 - ncore_;
      if (neleb < 0) throw runtime_error("Wrong states specified - there are not enough active electrons for the requested spin state.");
      if (nelea > norb_) throw runtime_error("Wrong states specified - there are not enough active orbitals for the requested spin state.");

      write_guess(nelea, neleb, nstate, offset);
      offset += nstate;
      if (nelea != neleb) {
        write_guess(neleb, nelea, nstate, offset);
        offset += nstate;
      }
    }
    pdebug.tick_print("guess generation");
  }*/    
  write_dice_input();
  //clean();
  dump_integrals();
  string dice_call = mpiprefix_ + " " + dice_dir_ + " input.dat >> output.dat" ;//+ "/RelDice";
  const char* command = dice_call.c_str();
  cout << command << endl;
  ifstream out_file("output.dat");
  if (!out_file) {
    system("rm output.dat");
    system("rm *.bkp");
    system("rm Dice*rdm*");
    system("rm shci.e");
  }
  system(command);
  //system("tail output.dat");
  //read energy from shci.e
  disk_energy();
}

// Currently not available after I removed const_denom
// This function is modified from ZHarrison::generate_guess
void ZSHCI::write_guess(ofstream& fs, const int nelea, const int neleb, const int nstate, const int offset) {
  shared_ptr<const Determinants> cdet = space_->finddet(nelea, neleb);
  int ndet = nstate*10;
  int oindex = offset;
  const bool spin_adapt = idata_->get<bool>("spin_adapt", true);
  while (oindex < offset+nstate) {
    cout << ndet << nelea << neleb;
    vector<pair<bitset<nbit__>, bitset<nbit__>>> bits = detseeds(ndet, nelea, neleb);
    // Spin adapt detseeds
    oindex = offset;
    vector<pair<bitset<nbit__>, bitset<nbit__>>> done;
    for (auto& it : bits) {
      bitset<nbit__> alpha = it.second;
      bitset<nbit__> beta = it.first;
      bitset<nbit__> open_bit = (alpha^beta);

      // This can happen if all possible determinants are checked without finding nstate acceptable ones.
      if (alpha.count() + beta.count() != nele_)
        throw logic_error("ZSHCI::write_guess produced an invalid determinant. Check the number of states being requested.");

      pair<bitset<nbit__>, bitset<nbit__>> config = spin_adapt ? make_pair(open_bit, alpha & beta) : it;
      if (find(done.begin(), done.end(), config) != done.end()) continue;
      done.push_back(config);

      // make sure that we have enough unpaired alpha
      const int unpairalpha = (alpha ^ (alpha & beta)).count();
      const int unpairbeta = (beta ^ (alpha & beta)).count();
      if (unpairalpha - unpairbeta < nelea - neleb) continue;
      pair<vector<tuple<int, int, int>>, double> adapt;
      if (spin_adapt) {
        adapt = space_->finddet(nelea, neleb)->spin_adapt(nelea-neleb, alpha, beta);
      } else {
        adapt.first = vector<tuple<int, int, int>>(1, make_tuple(space_->finddet(nelea, neleb)->lexical<1>(beta),
                                                                            space_->finddet(nelea, neleb)->lexical<0>(alpha), 1));
        adapt.second = 1.0;
      }
      cout << "print det" << endl;
      const double fac = adapt.second;
      vector<int> alphaIndex = bit_to_numbers(alpha);
      vector<int> betaIndex = bit_to_numbers(beta);
      for (auto& electron : alphaIndex)
        fs << 2*electron << " ";
      for (auto& electron : betaIndex)
        fs << 2*electron+1 << " ";
      fs << endl;
      ++oindex;
      if (oindex == offset+nstate) break;
    }
  }
}
std::shared_ptr<Kramers<2,ZRDM<1>>> ZSHCI::read_external_rdm1(const int ist, const int jst, const std::string& file) const {
  std::stringstream ss; ss << file << "_" << ist << "_" << jst << ".rdm1";
  std::ifstream in(ss.str(), std::ios::in | std::ios::binary);
  int norbs = 0;
  if (!in) {
    cout << "open error" << endl;
    exit(0);
  }
  in.read((char*) (&norbs), sizeof(int));
  if (norbs != 2*norb_)
    throw runtime_error("rdm given does not match with the current system");
  vector<complex<double>> data;
  data.resize(norbs*norbs); 
  in.read((char*) data.data(), norbs*norbs*sizeof(complex<double>));
  in.close();
  auto out = make_shared<Kramers<2,ZRDM<1>>>();
  auto tmp = make_shared<ZRDM<1>>(norb_);
  out->add(KTag<2>("00"), tmp->clone());
  out->add(KTag<2>("01"), tmp->clone());
  out->add(KTag<2>("10"), tmp->clone());
  out->add(KTag<2>("11"), tmp);
  for (int i=0; i<norbs; i++) {
    for (int j=0; j<norbs; j++) {
      const int jj = i/2;
      const int ii = j/2;
      const complex<double> dat = data[i*norbs+j];
      const int ti = j%2;
      const int tj = i%2;
      const KTag<2> tag1{ti, tj};
      out->at(tag1)->element(ii, jj) = dat;
      if (ist == jst) {
        const KTag<2> tag2{tj, ti};
        out->at(tag2)->element(jj, ii) = conj(dat);
      }
    }
  }
  return out;
}
std::shared_ptr<Kramers<4,ZRDM<2>>> ZSHCI::read_external_rdm2(const int ist, const int jst, const std::string& file) const {
  std::stringstream ss; ss << file << "_" << ist << "_" << jst << ".rdm2";
  std::ifstream in(ss.str(), std::ios::in | std::ios::binary);
  if (!in) {
    cout << "open error" << endl;
    exit(0);
  }
  int norbs = 0;

  in.read((char*) (&norbs), sizeof(int));
  cout << norbs << " " << norb_ << endl;
  if (norbs != 2*norb_)
    throw runtime_error("rdm given does not match with the current system");
  vector<complex<double>> data;
  int vec_size = pow(norbs, 4);
  data.resize(vec_size);
  in.read((char*) data.data(), vec_size*sizeof(complex<double>));
  in.close();
  
  auto out = make_shared<Kramers<4,ZRDM<2>>>();

  map<array<int,2>, double> elem;
  elem.emplace(array<int,2>{{0,1}}, 1.0); elem.emplace(array<int,2>{{1,0}},-1.0);
  #pragma omp for
  for (int ij=0; ij<norbs*norbs; ij++) {
    for (int kl=0; kl<norbs*norbs; kl++) {
      const int i=ij/norbs;
      const int j=ij%norbs;
      const int k=kl/norbs;
      const int l=kl%norbs;
      const complex<double> dat=data[ij*norbs*norbs+kl];

      map<int,pair<int,int>> mij{{0,{(i)/2,(i)%2}}, {1,{(j)/2,(j)%2}}};
      map<int,pair<int,int>> mkl{{0,{(k)/2,(k)%2}}, {1,{(l)/2,(l)%2}}};

      map<int, int> aij{{0,i},{1,j}};
      map<int, int> akl{{0,k},{1,l}};
      for (auto& eij : elem) {
        for (auto& ekl : elem) {
          //if (aij[eij.first[0]] < aij[eij.first[1]] || akl[ekl.first[0]] < akl[ekl.first[1]]) continue;
          if (mij[eij.first[0]].second > mij[eij.first[1]].second || mkl[ekl.first[0]].second > mkl[ekl.first[1]].second) continue;
          //if (abs(dat) < 1.0e-20) continue;
          const KTag<4> t{mij[eij.first[0]].second, mkl[ekl.first[0]].second, mij[eij.first[1]].second, mkl[ekl.first[1]].second};
          if (!out->exist(t))
            out->add(t, make_shared<ZRDM<2>>(norb_));
	  double sgn = 1.0;
	  if (aij[eij.first[0]] >= aij[eij.first[1]] && akl[ekl.first[0]] >= akl[ekl.first[1]]) sgn = 1.0;
	  else if (aij[eij.first[0]] < aij[eij.first[1]] && akl[ekl.first[0]] < akl[ekl.first[1]]) sgn = 1.0;
          else sgn = -1.0;
	  out->at(t)->element(mij[eij.first[0]].first, mkl[ekl.first[0]].first, mij[eij.first[1]].first, mkl[ekl.first[1]].first)
            = sgn * dat;//eij.second * ekl.second * dat;
          //cout << scientific <<  setprecision(7) << dat << setprecision(2) << " " << eij.second << " " << ekl.second << " " << i << " " << j << " " << k << " " << l << endl;
          if (ist == jst) {
            const KTag<4> t2{mkl[ekl.first[0]].second, mij[eij.first[0]].second, mkl[ekl.first[1]].second, mij[eij.first[1]].second};
            if (!out->exist(t2))
              out->add(t2, make_shared<ZRDM<2>>(norb_));
            out->at(t2)->element(mkl[ekl.first[0]].first, mij[eij.first[0]].first, mkl[ekl.first[1]].first, mij[eij.first[1]].first)
              = eij.second * ekl.second * conj(dat);
          }
        }
      }
    }
  }

  for (auto& i : elem) {
    for (auto& j : elem) {
      vector<int> perm(4);
      for (int k=0; k!=2; ++k) {
        perm[k*2]   = j.first[k]*2;
        perm[k*2+1] = i.first[k]*2+1;
      }
      out->emplace_perm(perm, j.second*i.second);
    }
  }
  return out;
}
void ZSHCI::write_dice_input() {

  vector<vector<int>> seed_dets;
  auto dets = shci_data_->get_child_optional("guess_determinant");
  if (dets) {
    for (auto it = dets->begin(); it != dets->end(); it++) {
      auto row = *it;
      seed_dets.push_back(row->get_vector<int>("", 0));
    }
  }
  else {
    for (int i=nele_/2; i>0; i--) {
      if (i>norb_ | nele_-i>norb_) continue;
      else if (i==nele_-i) {
        vector<int> idet;
        for (int j=0; j<i; j++) {
          idet.push_back(2*j);
          idet.push_back(2*j+1);
        }
        seed_dets.push_back(idet);
      }
      else {
        vector<int> alpha_det;
        for (int j=0; j<i; j++) {
          alpha_det.push_back(2*j);
        }
        for (int j=0; j<nele_-i; j++) {
          alpha_det.push_back(2*j+1);
        }
        vector<int> beta_det;
        for (int j=0; j<i; j++) {
          beta_det.push_back(2*j+1);
        }
        for (int j=0; j<nele_-i; j++) {
          beta_det.push_back(2*j);
        }
        seed_dets.push_back(alpha_det);
        seed_dets.push_back(beta_det);
      }
    }
  }
 
  vector<int> sweep_iter;
  vector<double> sweep_epsilon;

  if (shci_data_->get_child_optional("sweep_iter")) { 
    sweep_iter = shci_data_->get_vector<int>("sweep_iter");
    sweep_epsilon = shci_data_->get_vector<double>("sweep_epsilon");
  }
  else {
    sweep_iter = {0, 2};
    sweep_epsilon = {1e-2, 1e-5};
  }

  if (sweep_iter.size() != sweep_epsilon.size()) throw runtime_error("Size of sweep_iter and sweep_epsilon has to be the same");
  
  int maxiter = sweep_iter[sweep_iter.size()-1]+6;
  bool io = shci_data_->get<bool>("diskio", false);
  bool restart = shci_data_->get<bool>("restart", false);
  bool fullrestart = shci_data_->get<bool>("fullrestart", false);
  bool stochastic = shci_data_->get<bool>("stochasticPT", false);
  bool trev = shci_data_->get<bool>("Treversal", 0);
  int PTiter = shci_data_->get<int>("PTiter", 0);
  double davidsonTol = shci_data_->get<double>("davidson", 1e-7);
  double davidsonTolLoose = shci_data_->get<double>("davidsonLoose",  1e-5);
  double epsilon2 = shci_data_->get<double>("epsilon2", 1e-10);
  double epsilon2Large = shci_data_->get<double>("epsilon2Large", 1e-5);
  double targetError = shci_data_->get<double>("targetError", 1e-4);
  double dE = shci_data_->get<double>("dE", 1e-8);
  int nroot = nstate_;
  int sampleN = shci_data_->get<int>("sampleN", 100);

  ofstream fs("input.dat");
  fs << "#system" << endl;
  fs << "nocc " << nele_ << endl;

  for (auto & row : seed_dets) {
    for (auto & i : row) {
      fs << i << " ";
    }
    fs << endl;
  } 
  // generate the initial guess
  // partially copied from the generate_guess block in zharrison compute function
  // Temperarity remove this function due to the problem caused by removing const_denom.
  /*int offset = 0;
  for (int ispin = 0; ispin != states_.size(); ++ispin) {
    int nstate = 0;
    for (int i = ispin; i != states_.size(); ++i)
      nstate += states_[i];

    if (nstate == 0)
      continue;
    
    if ((geom_->nele()+ispin-charge_) % 2 == 1) {
      if (states_[ispin] == 0) {
        continue;
      } else {
        if ((geom_->nele()-charge_) % 2 == 0) throw runtime_error("Wrong states specified - only integer spins are allowed for even electron counts.");
        else throw runtime_error("Wrong states specified - only half-integer spins are allowed for odd electron counts.");
      }
    }

    const int nelea = (geom_->nele()+ispin-charge_)/2 - ncore_;
    const int neleb = (geom_->nele()-ispin-charge_)/2 - ncore_;
    if (neleb < 0) throw runtime_error("Wrong states specified - there are not enough active electrons for the requested spin state.");
    if (nelea > norb_) throw runtime_error("Wrong states specified - there are not enough active orbitals for the requested spin state.");
    cout << "write guess" << endl;
    write_guess(fs, nelea, neleb, nstate, offset);
    offset += nstate;
    if (nelea != neleb) {
      write_guess(fs, neleb, nelea, nstate, offset);
      offset += nstate;
    }
  }*/
  // Read guess determinant from input instead
  fs << "end" << endl;
  if (fullrestart && !first_iter_) {
    restart = true;
    sweep_iter = vector<int>(1,0);
    sweep_epsilon = vector<double>(1, 10.0);
  } 
  fs << "nroots " << nroot << endl;
  fs << "#variational" << endl;
  fs << "schedule" << endl;
  for (int i=0; i<sweep_iter.size(); i++)
    fs << sweep_iter[i] << " " << sweep_epsilon[i] << endl;
  fs << "end" << endl;
  fs << "davidsonTol " << davidsonTol << endl;
  fs << "davidsonTolLoose " << davidsonTolLoose << endl;
  fs << "dE " << dE << endl;
  fs << "maxiter " << maxiter << endl;
  if (restart && !first_iter_) {
    fs << "fullrestart" << endl;
  }  
  fs << "#pt" << endl;
  if (PTiter != 0 || tight_) {
    if(!stochastic)
      fs << "deterministic" << endl;
    fs << "epsilon2 " << epsilon2 << endl;
    fs << "epsilon2Large " << epsilon2Large << endl;
    fs << "targetError " << targetError << endl;
    fs << "sampleN " << sampleN << endl;
  }
  else
    fs << "nPTiter 0" << endl;

  fs << "#misc" << endl;
  if (!io)
    fs << "noio" << endl;
  if (!trev)
    fs << "Treversal " << trev << endl;
  if (!tight_) {
    fs << "DoSpinRDM" << endl;
    fs << "DoOneRDM" << endl;
  }
  fs << endl;
  fs.close();

  return;
}
void ZSHCI::compute_rdm12() {
  rdm1_.clear();
  rdm2_.clear();
  rdm1_.resize(nstate_);
  rdm2_.resize(nstate_);
  //read_dice_rdm();
  for (int istate = 0; istate != nstate_; ++istate) {
    // one body RDM
    rdm1_[istate] = read_external_rdm1(istate, istate, rdm_file_);

    // if off-diagonals are zero, generate a blank RDM for completeness
    if (!rdm1_[istate]->exist({1,0}))
      rdm1_[istate]->add({1,0}, rdm1_[istate]->at({0,0})->clone());

    // two body RDM
    rdm2_[istate] = read_external_rdm2(istate, istate, rdm_file_);

    // append permutation information
    //rdm2_[istate]->emplace_perm({{0,3,2,1}},-1);
    //rdm2_[istate]->emplace_perm({{2,3,0,1}}, 1);
    //rdm2_[istate]->emplace_perm({{2,1,0,3}},-1);
    
  }

  if (nstate_ > 1) {
    rdm1_av_ = make_shared<Kramers<2,ZRDM<1>>>();
    rdm2_av_ = make_shared<Kramers<4,ZRDM<2>>>();
    for (int istate = 0; istate != nstate_; ++istate) {
      for (auto& i : *rdm1_[istate])
        rdm1_av_->add(i.first, i.second);
      for (auto& i : *rdm2_[istate])
        rdm2_av_->add(i.first, i.second);
    }
    for (auto& i : *rdm1_av_) i.second->scale(1.0/nstate_);
    for (auto& i : *rdm2_av_) i.second->scale(1.0/nstate_);
    rdm2_av_->emplace_perm({{0,3,2,1}},-1);
    rdm2_av_->emplace_perm({{2,3,0,1}}, 1);
    rdm2_av_->emplace_perm({{2,1,0,3}},-1);
  } else {
    rdm1_av_ = rdm1_.front();
    rdm2_av_ = rdm2_.front();
  }
  // set expanded
  rdm1_av_expanded_ = expand_kramers<1,complex<double>>(rdm1_av_, norb_);
  rdm2_av_expanded_ = expand_kramers<2,complex<double>>(rdm2_av_, norb_);
  if (idata_->get<bool>("dump_rdms", false)) {
    const double rdm_thresh = idata_->get<double>("rdm_thresh_small",1.0e-5);
    rdm1_[0]->print();
    rdm2_[0]->print();
    rdm1_av_expanded_->print();
    rdm2_av_expanded_->print();
  }
  //clean();
  /*if (idata_->get<bool>("dump_rdms", true)) {
    auto tmp = make_shared<ZMatrix>(4*norb_*norb_, 4*norb_*norb_);
    copy_n(rdm2_av_expanded_->data(), tmp->size(), tmp->data());
    ofstream fs("DICE.rdm2");
    double upper = idata_->get<double>("rdm_thresh_large", 1.00);
    double thresh = idata_->get<double>("rdm_thresh_small", 0.05);
    for (int i = 0; i != norb_*2; ++i)
        for (int j = 0; j != norb_*2; ++j)
          for (int k = 0; k != norb_*2; ++k)
            for (int l = 0; l != norb_*2; ++l) {
              const complex<double> val = tmp->element(j+norb_*i, l+norb_*k);
              if (abs(val) > thresh && abs(val) < upper) {
                fs << setw(20) << val.real() << setw(20) << val.imag() << setw(4) << i+1 << setw(4)
                  << j+1 << setw(4) << k+1 << setw(4) << l+1 << endl;   // electron 1, electron 2
              }
            }
    fs.close();
    auto tmp2 = make_shared<ZMatrix>(2*norb_, 2*norb_);
    copy_n(rdm1_av_expanded_->data(), tmp2->size(), tmp2->data());
    ofstream fs2("DICE.rdm1");
    for (int i=0; i != norb_*2; ++i)
      for (int j=0; j != norb_*2; ++j) {
        const complex<double> val = tmp2->element(i,j);
        if (abs(val) > thresh && abs(val) < upper)
          fs2 << setw(20) << val << setw(4) << i+1 << setw(4) << j+1 << endl;
      }
    fs2.close();
  }*/
  return;
}
