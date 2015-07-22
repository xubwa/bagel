//
// BAGEL - Parallel electron correlation program.
// Filename: nevpt2_rdm.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifdef NEVPT2IMPL

template<>
void NEVPT2<double>::compute_rdm() {
  // rdm 1
  {
    auto tmp = ref_->rdm1(istate_)->rdm1_mat(/*nclosed_*/0, false);
    tmp->localize();
    rdm1_ = tmp;
  }
  // rdm 2
  {
    auto tmp = make_shared<MatType>(nact_*nact_, nact_*nact_, true);
    shared_ptr<const RDM<2>> r2 = ref_->rdm2(istate_);
    sort_indices<0,2,1,3,0,1,1,1>(r2->data(), tmp->data(), nact_, nact_, nact_, nact_);
    rdm2_ = tmp;
  }
  // rdm 3 and 4
  {
    auto tmp3 = make_shared<MatType>(nact_*nact_*nact_, nact_*nact_*nact_, true);
    auto tmp4 = make_shared<MatType>(nact_*nact_*nact_*nact_, nact_*nact_*nact_*nact_, true);
    shared_ptr<const RDM<3>> r3;
    shared_ptr<const RDM<4>> r4;
    tie(r3, r4) = ref_->rdm34(istate_, istate_);
    sort_indices<0,2,4,  1,3,5,  0,1,1,1>(r3->data(), tmp3->data(), nact_, nact_, nact_, nact_, nact_, nact_);
    sort_indices<0,2,4,6,1,3,5,7,0,1,1,1>(r4->data(), tmp4->data(), nact_, nact_, nact_, nact_, nact_, nact_, nact_, nact_);
    rdm3_ = tmp3;
    rdm4_ = tmp4;
  }
}

template<>
void NEVPT2<complex<double>>::compute_rdm() {
  auto ref = dynamic_pointer_cast<const RelReference>(ref_);
  // rdm 1
  {
    auto tmp = make_shared<MatType>(nact_, nact_, true);
    auto rdm1k = ref->rdm1(istate_, istate_);
    auto r1 = expand_kramers(rdm1k, nact_/2);
    copy_n(r1->data(), nact_*nact_, tmp->data());
    rdm1_ = tmp;
  }
  // rdm 2
  {
    auto tmp = make_shared<MatType>(nact_*nact_, nact_*nact_, true);
    auto rdm2k = ref->rdm2(istate_, istate_);
    auto r2 = expand_kramers(rdm2k, nact_/2);
    sort_indices<0,2,1,3,0,1,1,1>(r2->data(), tmp->data(), nact_, nact_, nact_, nact_);
    rdm2_ = tmp;
  }
  // rdm 3
  {
    auto tmp3 = make_shared<MatType>(nact_*nact_*nact_, nact_*nact_*nact_, true);
    auto rdm3k = ref->rdm3(istate_, istate_);
    auto r3 = expand_kramers(rdm3k, nact_/2);
    sort_indices<0,2,4,1,3,5,0,1,1,1>(r3->data(), tmp3->data(), nact_, nact_, nact_, nact_, nact_, nact_);
    rdm3_ = tmp3;
  }
  // TODO rdm 4 is too large to do this way - implement direct computation in ARDM3 later
  // rdm 4
  {
    auto tmp4 = make_shared<MatType>(nact_*nact_*nact_*nact_, nact_*nact_*nact_*nact_, true);
    auto rdm4k = ref->rdm4(istate_, istate_);
    auto r4 = expand_kramers(rdm4k, nact_/2);
    sort_indices<0,2,4,6,1,3,5,7,0,1,1,1>(r4->data(), tmp4->data(), nact_, nact_, nact_, nact_, nact_, nact_, nact_, nact_);
    rdm4_ = tmp4;
  }
}


template<typename DataType>
void NEVPT2<DataType>::compute_asrdm() {
  assert(rdm1_ && rdm2_ && rdm3_ && rdm4_);
  auto id2 = [this](                          const int k, const int l) { return         (        (k+nact_*l)); };
  auto id3 = [this](             const int j, const int k, const int l) { return         (j+nact_*(k+nact_*l)); };
  auto id4 = [this](const int i, const int j, const int k, const int l) { return i+nact_*(j+nact_*(k+nact_*l)); };

  const double fac2 = is_same<DataType,double>::value ? 2.0 : 1.0;

  shared_ptr<MatType> srdm2 = rdm2_->clone(); // S(a,b,c,d) = <0|a+p bp cq d+q|0>
  for (int i = 0; i != nact_; ++i)
    for (int j = 0; j != nact_; ++j)
      for (int k = 0; k != nact_; ++k)
        for (int l = 0; l != nact_; ++l)
          srdm2->element(l+nact_*k,j+nact_*i) = -rdm2_->element(l+nact_*i,k+nact_*j) + (i == j ? fac2*rdm1_->element(l,k) : 0.0) - (i == k ? rdm1_->element(l,j) : 0.0);
  // <a+ a b+ b> and <a+ a b+ b c+ c>
  shared_ptr<MatType> ardm2 = rdm2_->clone();
  for (int i = 0; i != nact_; ++i)
    for (int j = 0; j != nact_; ++j)
      for (int k = 0; k != nact_; ++k) {
        for (int l = 0; l != nact_; ++l)
          ardm2->element(l+nact_*k,j+nact_*i) += rdm2_->element(l+nact_*j,k+nact_*i);
        ardm2->element(k+nact_*j,j+nact_*i) += rdm1_->element(k,i);
      }
  shared_ptr<MatType> ardm3 = rdm3_->clone();
  shared_ptr<MatType> srdm3 = rdm3_->clone(); // <a+ a b b+ c+ c>
  for (int i = 0; i != nact_; ++i)
    for (int j = 0; j != nact_; ++j)
      for (int k = 0; k != nact_; ++k)
        for (int l = 0; l != nact_; ++l)
          for (int m = 0; m != nact_; ++m) {
            for (int n = 0; n != nact_; ++n)
              ardm3->element(id3(n,m,l),id3(k,j,i)) += rdm3_->element(id3(n,l,j),id3(m,k,i));
            ardm3->element(id3(m,l,l),id3(k,j,i)) += ardm2->element(m+nact_*k,j+nact_*i);
            ardm3->element(id3(m,l,k),id3(j,j,i)) += rdm2_->element(m+nact_*k,l+nact_*i);
            ardm3->element(id3(m,l,k),id3(j,l,i)) += rdm2_->element(m+nact_*k,i+nact_*j);

            srdm3->element(id3(m,l,k),id3(k,j,i)) += fac2*ardm2->element(id2(m,l),id2(j,i));
          }
  sort_indices<0,2,1,3,1,1,-1,1>(ardm3->data(), srdm3->data(), nact_*nact_, nact_, nact_, nact_*nact_);
  shared_ptr<MatType> ardm4 = rdm4_->clone();
  for (int h = 0; h != nact_; ++h)
    for (int g = 0; g != nact_; ++g)
      for (int f = 0; f != nact_; ++f)
        for (int e = 0; e != nact_; ++e)
          for (int d = 0; d != nact_; ++d)
            for (int c = 0; c != nact_; ++c)
              for (int b = 0; b != nact_; ++b)
                for (int a = 0; a != nact_; ++a) {
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (b == c ? 1.0 : 0.0) * ardm3->element(id3(a,d,e),id3(f,g,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) -= (d == e && b == c ? 1.0 : 0.0) * ardm2->element(id2(a,f),id2(g,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (d == e ? 1.0 : 0.0) * ardm3->element(id3(a,b,c),id3(f,g,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) -= (b == e && c == f ? 1.0 : 0.0) * ardm2->element(id2(a,d),id2(g,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (b == e ? 1.0 : 0.0) * ardm3->element(id3(a,f,c),id3(d,g,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (f == g ? 1.0 : 0.0) * rdm3_->element(id3(a,c,e),id3(b,d,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (d == g ? 1.0 : 0.0) * rdm3_->element(id3(a,c,e),id3(b,h,f));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (b == g ? 1.0 : 0.0) * rdm3_->element(id3(a,c,e),id3(h,d,f));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += rdm4_->element(id4(a,c,e,g),id4(b,d,f,h));
                }
  ardm2_ = ardm2;
  ardm3_ = ardm3;
  ardm4_ = ardm4;
  srdm2_ = srdm2;
  srdm3_ = srdm3;
}


template<typename DataType>
void NEVPT2<DataType>::compute_hrdm() {
  assert(rdm1_ && rdm2_ && rdm3_ && srdm2_);

  auto id3 = [this](const int j, const int k, const int l) { return j+nact_*(k+nact_*l); };

  const double fac2 = is_same<DataType,double>::value ? 2.0 : 1.0;

  shared_ptr<MatType> unit = rdm1_->clone(); unit->unit();
  shared_ptr<MatType> hrdm2 = rdm2_->copy();
  shared_ptr<const MatType> hrdm1 = make_shared<MatType>(*unit*fac2 - *rdm1_);

  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nact_; ++j) {
      for (int k = 0; k != nact_; ++k) {
        hrdm2->element(j+nact_*k, i+nact_*k) += fac2 * hrdm1->element(j,i);
        hrdm2->element(k+nact_*j, i+nact_*k) -= hrdm1->element(j,i);
        hrdm2->element(j+nact_*k, k+nact_*i) += rdm1_->element(j,i);
        hrdm2->element(k+nact_*j, k+nact_*i) -= fac2 * rdm1_->element(j,i);
      }
    }
  }
  auto hrdm3 = make_shared<MatType>(*rdm3_ * (-1.0));
  for (int i = 0; i != nact_; ++i)
    for (int j = 0; j != nact_; ++j)
      for (int k = 0; k != nact_; ++k)
        for (int l = 0; l != nact_; ++l)
          for (int m = 0; m != nact_; ++m) {
            hrdm3->element(id3(l,k,m),id3(j,i,m)) += fac2*hrdm2->element(l+nact_*k,j+nact_*i);
            hrdm3->element(id3(l,m,k),id3(j,i,m)) -=      hrdm2->element(l+nact_*k,j+nact_*i);
            hrdm3->element(id3(m,l,k),id3(j,i,m)) -=      hrdm2->element(l+nact_*k,i+nact_*j);
            hrdm3->element(id3(l,k,m),id3(j,m,i)) +=      detail::conj(srdm2_->element(i+nact_*k,l+nact_*j));
            hrdm3->element(id3(l,m,k),id3(j,m,i)) -= fac2*detail::conj(srdm2_->element(i+nact_*k,l+nact_*j));
            hrdm3->element(id3(m,l,k),id3(j,m,i)) +=      detail::conj(srdm2_->element(i+nact_*k,l+nact_*j));
            hrdm3->element(id3(l,k,m),id3(m,j,i)) -=      rdm2_->element(l+nact_*k,i+nact_*j);
            hrdm3->element(id3(l,m,k),id3(m,j,i)) -=      rdm2_->element(k+nact_*l,i+nact_*j);
            hrdm3->element(id3(m,l,k),id3(m,j,i)) += fac2*rdm2_->element(k+nact_*l,i+nact_*j);
          }
  hrdm1_ = hrdm1;
  assert(hrdm1_->is_hermitian());
  hrdm2_ = hrdm2;
  assert(hrdm2_->is_hermitian());
  hrdm3_ = hrdm3;
  assert(hrdm3_->is_hermitian());
}

#endif