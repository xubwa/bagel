//
// BAGEL - Parallel electron correlation program.
// Filename: stevensop.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynolds2018@u.northwestern.edu>
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

// A function to compute the coefficients of Extended Stevens Operators, for the pseudospin Hamiltonian
// Notation follows I. D. Ryabov, Appl. Magn. Reson. (2009) 35, 481-494.
// Some equations also come from I. D. Ryabov, J. Magn. Reson. (1999) 140, 141-145.

#include <src/prop/pseudospin/pseudospin.h>
#include <src/util/math/factorial.h>
#include <src/util/math/comb.h>
#include <src/util/math/matop.h>


using namespace std;
using namespace bagel;

// Some helper functions only used here
namespace {

const Factorial fact;
const Comb comb;

// Even k:  alpha = 1.0
// Odd k:   alpha = 1.0 or 0.5 for even or odd q, respectively
double compute_alpha(const int k, const int q) {
  double out;
  if (k % 2 == 0)
    out = 1.0;
  else
    out = (q % 2 == 0) ? 1.0 : 0.5;
  return out;
}

// Compute Nkq for q = 0
double compute_Nk0(const int k) {
  assert (k >= 0);
  const double twok = std::pow(2.0,k);
  const double Nkk_denomenator = twok * fact(k);
  const double out = 1.0 / Nkk_denomenator;
  return out;
}

// a(k, q; m, i) in Ryabov's notation
// Uses a downward recurrence relation:  the input vector akq1m contains the output for k, q+1, and all values of m, i
// Output corresponds to a specific k, q; m but all values of i
vector<long long> compute_akqmi(const int k, const int q, const int m, const int nspin, const vector<vector<long long>> akq1mi) {
  assert (k >= 1 && q >= 0 && m >= 0 && q < k && m <= k - q);

  const int nval = (k - q - m) / 2 + 1;
  vector<long long> out(nval);

  for (int i=0; i!=nval; ++i) {

    out[i] = m > 0 ? (2 * q + m + 1) * akq1mi[m - 1][i] : 0;

    if (m <= k - q - 1 && i <= (k - q - m - 1) / 2)
      out[i] += (q * (q + 1) - m * (m + 1) / 2) * akq1mi[m][i];

    for (int n = 1; n <= k - q - m - 1; ++n) {
      const long long sign = (n % 2 == 0) ? 1 : -1;
      const long long coeff1 = comb(m + n, m);
      const long long coeff2 = m > 0 ? comb(m + n, m - 1) : 0;
      const long long coeff3 = m > 1 ? comb(m + n, m - 2) : 0;

      if (i <= (k - q - m - n + 1) / 2 && i > 0)
        out[i] += sign * coeff1 * akq1mi[m + n][i-1];

      if (i <= (k - q - m - n - 1) / 2)
        out[i] -= sign * (coeff2 + coeff3) * akq1mi[m + n][i];
    }
  }
  return out;
}

long long greatest_common_factor(const long long a, const long long b) {
  return b == 0 ? a : greatest_common_factor(b, a % b);
}

// Fkq is the greatest common factor of all a(k, q; m, i) with the same k, q
long long compute_Fkq(const vector<vector<long long>> input) {
  long long out = std::abs(input[0][0]);
  for(int m = 0; m != input.size(); ++m)
    for(int i = 0; i != input[m].size(); ++i)
      out = greatest_common_factor(out, std::abs(input[m][i]));
  return out;
}

} // end of anonymous namespace


// And now the driver
vector<Spin_Operator> Pseudospin::build_extended_stevens_operators(const vector<int> ranks) const {
  cout << "    Computing extended stevens operators for S = " << nspin_/2 << (nspin_ % 2 == 0 ? "" : " 1/2") << endl;
  const double ss1 = nspin_ * (nspin_ + 2.0) / 4.0; // S(S+1)

  vector<Spin_Operator> stevensop = {};
  cout << fixed << setprecision(6);

  for (auto& k : ranks) {

    // Requires factorial of k
    // Above k = 12, coefficients become so large that long long fails to capture them.
    // Potentially we could store a(k, q; m) / Fkq or something, but for now this is fine.
    const int kmax = std::min(fact.max(), 13);
    if (k >= kmax)
      throw runtime_error("Sorry, numerical issues currently limit us to Stevens operators of order " + to_string(kmax - 1) + " and lower");
    if (k < 0)
      throw runtime_error("Ranks of Extended Stevens Operators must be whole numbers.");

    vector<double> alpha(k + 1);
    vector<double> Nkq(k + 1);  // positive q
    vector<double> Nk_q(k + 1); // negative q

    // a(k, q; m) in Ryabov's notation
    // access is akqm[q][m]
    vector<vector<double>> akqm(k + 1);

    // same but for negative q
    vector<vector<double>> ak_qm(k + 1);

    // The coefficients that go into computation of akqm
    vector<vector<vector<long long>>> akqmi(k + 1);

    alpha[0] = 1.0;
    Nkq[0] = compute_Nk0(k);

    for (int q = 1; q <= k; ++q) {
      alpha[q] = compute_alpha(k, q);
      Nkq[q] = -1.0 * Nkq[q-1] * std::sqrt((k + q) * (k - q + 1.0));
    }

    for (int q = k; q >= 0; --q) {

      vector<double> akqm_current(k - q + 1);  // positive q
      vector<double> ak_qm_current(k - q + 1); // negative q

      vector<vector<long long>> akqmi_current(k - q + 1);

      Nk_q[q] = Nkq[q] * ((k % 2 == 0) ? 1.0 : -1.0);

      for (int m = 0; m <= k - q; ++m) {
        if (q == k)
          akqmi_current[m] = vector<long long>(1, 1);
        else
          akqmi_current[m] = compute_akqmi(k, q, m, nspin_, akqmi[q + 1]);

        akqm_current[m] = 0.0;
        for (int i=0; i <= (k-q-m)/2; ++i)
          akqm_current[m] += std::pow(ss1, i) * akqmi_current[m][i];

        ak_qm_current[m] = akqm_current[m] * (m % 2 == 0 ? 1.0 : -1.0);
      }

      const long long Fkq = compute_Fkq(akqmi_current);

      // Debug printout
      //const double Nkk = Nkq[0]*std::sqrt(fact(2*k))*(k % 2 == 0 ? 1.0 : -1.0);
      //for (int m = 0; m <= k - q; ++m)
        //cout << "       k = " << setw(2) << k << ", q = " << setw(2) << q << ", m = " << setw(2) << m << ", alpha = " << setw(6) << setprecision(2) << alpha[q] << ", Nkk = " << setw(10) << setprecision(6) << Nkq[0]*std::sqrt(fact(2*k))*(k % 2 == 0 ? 1.0 : -1.0) << ", Nkq = " << setw(10) << Nkq[q] << ", Nk_q = " << setw(10) << Nk_q[q] << ", Nkq/Nkk = " << setw(10) << Nkq[q] / Nkk << "    akqm = " << setw(14) << setprecision(2) << akqm_current[m] << ", Fkq = " << setw(6) << Fkq << ", a(k, q; m)/Fkq = " << setw(10) << akqm_current[m] / Fkq << endl;

      akqm[q] = akqm_current;
      akqmi[q] = akqmi_current;
      ak_qm[q] = ak_qm_current;

      shared_ptr<ZMatrix> Tkq = spin_plus()->clone();
      shared_ptr<ZMatrix> Tk_q = spin_plus()->clone();
      shared_ptr<ZMatrix> Tk_q_check = spin_plus()->clone();
      const double sign1 = (q % 2 == 0) ? 1.0 : -1.0;

      for (int m = 0; m <= k - q; ++m) {
        // Need spin-z matrix to power of m and spin-plus to power of q
        shared_ptr<ZMatrix> Szm = spin_xyz(2)->clone();
        shared_ptr<ZMatrix> Spq = spin_xyz(2)->clone();
        Szm->unit();
        Spq->unit();
        for (int i = 0; i != m; ++i)
          *Szm *= *spin_xyz(2);
        for (int i = 0; i != q; ++i)
          *Spq *= *spin_plus();

        const double sign2 = ((k - m) % 2 == 0) ? 1.0 : -1.0;
        const double coeff = sign1 * sign2 * Nkq[q] * akqm[q][m];
        *Tkq += (coeff * *Szm * *Spq);

        /******/ // Debug code
        shared_ptr<ZMatrix> Smq = spin_xyz(2)->clone();
        Smq->unit();
        for (int i = 0; i != q; ++i)
          *Smq *= *spin_minus();
        const double coeff2 = sign1 * Nkq[q] * akqm[q][m];
        *Tk_q_check += (coeff2 * *Szm * *Smq);
        /******/
      }

      *Tk_q = sign1 * *Tkq->transpose_conjg();

      assert(((*Tk_q - *Tk_q_check).rms() < Tk_q->rms() * 1.0e-8) || (Tk_q->rms() < 1.0e-8));

      const double ckq = alpha[q] / (Nkq[q] * Fkq);
      auto Okq = make_shared<const ZMatrix>(ckq * *Tkq);
      auto Ok_q = make_shared<const ZMatrix>(sign1 * ckq * *Tk_q);

      // assert that the +q and -q algorithms give the same final result for q = 0
      assert(q != 0 || (*Okq - *Ok_q).rms() < 1.0e-10);

      auto Ocos_kq = make_shared<const ZMatrix>(sign1 * 0.5 * (*Okq + *Okq->transpose_conjg()));
      auto Osin_kq = make_shared<const ZMatrix>(sign1 * complex<double>(0.0, -0.5) * (*Okq - *Okq->transpose_conjg()));

      const string label_pq = "B_" + to_string(k) + "^" + to_string(q);
      const string label_mq = "B_" + to_string(k) + "^-" + to_string(q);
      stevensop.push_back(Spin_Operator(Ocos_kq, label_pq));
      if (q != 0)
        stevensop.push_back(Spin_Operator(Osin_kq, label_mq));

      const string spinstring = to_string(nspin_ / 2) + (nspin_ % 2 == 0 ? "" : " 1/2");
      if (Ocos_kq->rms() > 1.0e-10)
        Ocos_kq->print("Stevens operator, S = " + spinstring + ", k = " + to_string(k) + " , q = " + to_string(q), 12);
      if (Osin_kq->rms() > 1.0e-10)
        Osin_kq->print("Stevens operator, S = " + spinstring + ", k = " + to_string(k) + " , q = " + to_string(-q), 12);

    }
    cout << endl;
  }
  return stevensop;
}


