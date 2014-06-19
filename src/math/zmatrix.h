//
// BAGEL - Parallel electron correlation program.
// Filename: zmatrix.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Matthew Kelley <matthewkelley2017@u.northwestern.edu>
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


#ifndef __SRC_UTIL_ZMATRIX_H
#define __SRC_UTIL_ZMATRIX_H

#include <complex>
#include <src/math/matrix.h>

namespace bagel {

#ifdef HAVE_SCALAPACK
class DistZMatrix;
#else
class ZMatrix;
using DistZMatrix = ZMatrix;
#endif

class ZMatrix : public Matrix_base<std::complex<double>>, public std::enable_shared_from_this<ZMatrix> {
  public:
#ifdef HAVE_SCALAPACK
    ZMatrix(const int n, const int m, const bool localized = false);
#else
    ZMatrix(const int n, const int m, const bool localized = true);
#endif
    ZMatrix(const ZMatrix&);
    ZMatrix(const ZMatView& o);
    ZMatrix(ZMatrix&&);
    ZMatrix(const Matrix& real, const Matrix& imag);
    ZMatrix(const Matrix& real, const std::complex<double> factor);
    ZMatrix() { }
    virtual ~ZMatrix() { }

    void antisymmetrize();
    void hermite();
    std::shared_ptr<ZMatrix> cut(const int nstart, const int nend) const { return get_submatrix(nstart, 0, nend-nstart, mdim()); }
    std::shared_ptr<ZMatrix> slice_copy(const int mstart, const int mend) const { return get_submatrix(0, mstart, ndim(), mend-mstart); }
    std::shared_ptr<ZMatrix> resize(const int n, const int m) const { return this->resize_impl<ZMatrix>(n, m); }
    std::shared_ptr<ZMatrix> merge(const std::shared_ptr<const ZMatrix> o) const { return this->merge_impl<ZMatrix>(o); }

    std::shared_ptr<ZMatView> slice(const int mstart, const int mend) const {
      auto low = {0, mstart};
      auto up  = {ndim(), mend};
      return std::make_shared<ZMatView>(this->range().slice(low, up), this->storage(), localized_);
    }

    // diagonalize this matrix (overwritten by a coefficient matrix)
    virtual void diagonalize(double* vec);

    std::shared_ptr<ZMatrix> diagonalize_blocks(double* eig, std::vector<int> blocks) { return diagonalize_blocks_impl<ZMatrix>(eig, blocks); }

    std::tuple<std::shared_ptr<ZMatrix>, std::shared_ptr<ZMatrix>> svd(double* sing = nullptr);
    // compute S^-1. Assumes positive definite matrix
    void inverse();
    // compute S^-1/2. If an eigenvalue of S is smaller than thresh, the root will be discarded.
    bool inverse_half(const double thresh = 1.0e-8);
    std::shared_ptr<ZMatrix> tildex(const double thresh = 1.0e-8) const;

    using Matrix_base<std::complex<double>>::copy_block;
    using Matrix_base<std::complex<double>>::add_block;

    void copy_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const double* data);
    void copy_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const std::unique_ptr<double[]> o);
    void copy_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const std::shared_ptr<const Matrix> o);
    void copy_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const std::shared_ptr<const MatView> o);
    void copy_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const Matrix& o);

    void add_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const double* data);
    void add_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const std::unique_ptr<double[]> o);
    void add_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const std::shared_ptr<const Matrix> o);
    void add_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const std::shared_ptr<const MatView> o);
    void add_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const Matrix& o);

    std::shared_ptr<Matrix> get_real_part() const;
    std::shared_ptr<Matrix> get_imag_part() const;

    std::shared_ptr<ZMatrix> get_conjg() const;
    void fill_upper_conjg();

    std::shared_ptr<ZMatrix> get_submatrix(const int nstart, const int mstart, const int ndim, const int mdim) const {
      return this->get_submatrix_impl<ZMatrix>(nstart, mstart, ndim, mdim);
    }

    ZMatrix& operator=(const ZMatrix& o) { Matrix_base<std::complex<double>>::operator=(o); return *this; }
    ZMatrix& operator=(ZMatrix&& o)      { Matrix_base<std::complex<double>>::operator=(o); return *this; }

    ZMatrix& operator/=(const ZMatrix&);
    ZMatrix operator/(const ZMatrix&) const;

    std::shared_ptr<ZMatrix> clone() const { return std::make_shared<ZMatrix>(ndim(), mdim()); }
    std::shared_ptr<ZMatrix> copy() const { return std::make_shared<ZMatrix>(*this); }

    // returns exp(*this)
    std::shared_ptr<ZMatrix> exp(const int deg = 6) const;
    // returns log(*this)
    std::shared_ptr<ZMatrix> log(const int deg = 6) const;
    // returns transpose(*this)
    std::shared_ptr<ZMatrix> transpose() const;
    // returns hermite-conjugate(*this)
    std::shared_ptr<ZMatrix> transpose_conjg() const;

    using Matrix_base<std::complex<double>>::ax_plus_y;
    using Matrix_base<std::complex<double>>::dot_product;
    void ax_plus_y(const std::complex<double> a, const ZMatrix& o) { this->ax_plus_y_impl(a, o); }
    std::complex<double> dot_product(const ZMatrix& o) const { return this->dot_product_impl(o); }

    double orthog(const std::list<std::shared_ptr<const ZMatrix>> o) { return this->orthog_impl(o); }

    void add_diag(const std::complex<double> a, const int i, const int j) {
      assert(ndim() == mdim());
      for (int ii = i; ii != j; ++ii) element(ii,ii) += a;
    }
    void add_diag(const std::complex<double> a) { add_diag(a,0,ndim()); }
    // returns diagonal elements
    std::unique_ptr<std::complex<double>[]> diag() const;

    // purify a (near unitary) matrix to be unitary
    void purify_unitary();
    void purify_idempotent(const ZMatrix& s);
    void purify_redrotation(const int nclosed, const int nact, const int nvirt);

    std::shared_ptr<ZMatrix> solve(std::shared_ptr<const ZMatrix> A, const int n) const;

    void print(const std::string in = "", const int size = 10) const;

#ifdef HAVE_SCALAPACK
    std::shared_ptr<DistZMatrix> distmatrix() const;
    ZMatrix(const DistZMatrix&);
#else
    std::shared_ptr<const ZMatrix> distmatrix() const;
    std::shared_ptr<const ZMatrix> matrix() const { return shared_from_this(); }
#endif
    std::shared_ptr<const ZMatrix> form_density_rhf(const int n, const int off = 0, const std::complex<double> scale = 1.0) const;
};


#ifdef HAVE_SCALAPACK
// Not to be confused with Matrix. DistMatrix is distributed and only supported when SCALAPACK is turned on. Limited functionality
class DistZMatrix : public DistMatrix_base<std::complex<double>> {
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      boost::serialization::split_member(ar, *this, file_version);
    }
    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      std::shared_ptr<ZMatrix> mat = matrix();
      ar << mat;
    }
    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      std::shared_ptr<ZMatrix> mat;
      ar >> mat;
      DistZMatrix tmp(*mat);
      ndim_ = tmp.ndim_;
      mdim_ = tmp.mdim_;
      desc_ = tmp.desc_;
      localsize_ = tmp.localsize_;
      local_ = std::unique_ptr<std::complex<double>[]>(new std::complex<double>[tmp.size()]);
      *this = tmp;
    }

  public:
    DistZMatrix() { }
    DistZMatrix(const int n, const int m);
    DistZMatrix(const DistZMatrix&);
    DistZMatrix(DistZMatrix&&);
    DistZMatrix(const ZMatrix&);

    void diagonalize(double* vec) override;

    DistZMatrix operator*(const DistZMatrix&) const;
    DistZMatrix& operator*=(const DistZMatrix&);
    DistZMatrix operator%(const DistZMatrix&) const; // caution
    DistZMatrix operator^(const DistZMatrix&) const; // caution
    DistZMatrix operator+(const DistZMatrix& o) const { DistZMatrix out(*this); out.ax_plus_y(1.0, o); return out; }
    DistZMatrix operator-(const DistZMatrix& o) const { DistZMatrix out(*this); out.ax_plus_y(-1.0, o); return out; }
    DistZMatrix& operator+=(const DistZMatrix& o) { ax_plus_y(1.0, o); return *this; }
    DistZMatrix& operator-=(const DistZMatrix& o) { ax_plus_y(-1.0, o); return *this; }
    DistZMatrix& operator=(const DistZMatrix& o);
    DistZMatrix& operator=(DistZMatrix&& o);

    std::shared_ptr<DistZMatrix> clone() const { return std::make_shared<DistZMatrix>(ndim_, mdim_); }

    using DistMatrix_base<std::complex<double>>::scale;
    using DistMatrix_base<std::complex<double>>::ax_plus_y;
    using DistMatrix_base<std::complex<double>>::dot_product;

    void ax_plus_y(const std::complex<double> a, const DistZMatrix& o) { this->ax_plus_y_impl(a, o); }
    std::complex<double> dot_product(const DistZMatrix& o) const { return this->dot_product_impl(o); }

    std::shared_ptr<ZMatrix> matrix() const;

    std::shared_ptr<const DistZMatrix> form_density_rhf(const int n, const int off = 0, const std::complex<double> scale = 1.0) const;
};
#endif

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::ZMatrix)
#ifdef HAVE_SCALAPACK
BOOST_CLASS_EXPORT_KEY(bagel::DistZMatrix)
#endif

namespace bagel {
  template <class T>
  struct base_of<T, typename std::enable_if<std::is_base_of<ZMatrix, T>::value>::type> {
    typedef ZMatrix type;
  };
}

#endif
