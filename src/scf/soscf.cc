//
// BAGEL - Parallel electron correlation program.
// Filename: soscf.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <src/scf/soscf.h>
#include <src/scf/sofock.h>
#include <src/math/diis.h>
#include <src/scf/atomicdensities.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(SOSCF)

SOSCF::SOSCF(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
 : SCF_base(idata, geom, re) {
  cout << indent << "*** Two-component ECP-SCF ***" << endl << endl;
  soeig_ = unique_ptr<double[]> (new double[geom_->nbasis() * 2]);
  sohcore_base_ = make_shared<const SOHcore_base>(geom);

  sohcore_ = make_shared<SOHcore>(geom_, sohcore_base_);

}

void SOSCF::initial_guess() {
  sooverlap_ = sooverlap();
  sotildex_ = sotildex();

  shared_ptr<const ZMatrix> sofock = sohcore_;
  shared_ptr<ZMatrix> intermediate = make_shared<ZMatrix>(*sotildex_ % *sofock * *sotildex_);
  intermediate->diagonalize(soeig_.get());
  socoeff_ = make_shared<ZMatrix>(*sotildex_ * *intermediate);
  aodensity_ = socoeff_->form_density_rhf(nocc_ * 2);
}

void SOSCF::compute() {
  Timer scftime;
  initial_guess();

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl << endl;
  scftime.tick_print("SOSCF startup");
  cout << endl;
  cout << indent << "=== SOSCF iteration (" + geom_->basisfile() + ") ===" << endl << indent << endl;

  DIIS<DistZMatrix, ZMatrix> diis(diis_size_);

  for (int iter = 0; iter != max_iter_; ++iter) {
    shared_ptr<const ZMatrix> sofock = make_shared<const SOFock> (geom_, sohcore_, socoeff_->slice(0, nocc_ * 2));
    const complex<double> energy = 0.25 * ((*sohcore_ + *sofock) * *aodensity_).trace() + geom_->nuclear_repulsion();
    assert(energy.imag() < 1e-12);
    energy_ = energy.real();
    auto error_vector = make_shared<const ZMatrix>(*sofock * *aodensity_ * *sooverlap_ - *sooverlap_ * *aodensity_ * *sofock);
    auto real_error_vector = error_vector->get_real_part();
    const double error = real_error_vector->rms();

    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_ << "   "
                                      << setw(17) << error << setw(15) << setprecision(2) << scftime.tick() << endl;
    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * SOSCF iteration converged." << endl << endl;
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SOSCF." << endl << endl;
      break;
    }

    if (iter >= diis_start_) sofock = diis.extrapolate(make_pair(sofock, error_vector));

    shared_ptr<ZMatrix> intermediate = make_shared<ZMatrix>(*sotildex_ % *sofock * *sotildex_);
    intermediate->diagonalize(soeig_.get());
    socoeff_ = make_shared<ZMatrix>(*sotildex_ * *intermediate);
    aodensity_ = socoeff_->form_density_rhf(nocc_ * 2);
  }
}

shared_ptr<const ZMatrix> SOSCF::sotildex() {
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(2 * tildex_->ndim(), 2 * tildex_->mdim());
  out->zero();
  out->add_real_block(complex<double>(1.0, 0.0), 0, 0, tildex_->ndim(), tildex_->mdim(), tildex_);
  out->add_real_block(complex<double>(1.0, 0.0), tildex_->ndim(), tildex_->mdim(), tildex_->ndim(), tildex_->mdim(), tildex_);
  return out;
}

shared_ptr<const ZMatrix> SOSCF::sooverlap() {
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(2 * overlap_->ndim(), 2 * overlap_->mdim());
  out->zero();
  out->add_real_block(complex<double>(1.0, 0.0), 0, 0, overlap_->ndim(), overlap_->mdim(), overlap_);
  out->add_real_block(complex<double>(1.0, 0.0), overlap_->ndim(), overlap_->mdim(), overlap_->ndim(), overlap_->mdim(), overlap_);
  return out;
}

