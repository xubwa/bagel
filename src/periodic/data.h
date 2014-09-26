//
// BAGEL - Parallel electron correlation program.
// Filename: data.h
// Copyright (C) 2014 Toru Shiozaki
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


#ifndef __SRC_PERIODIC_DATA_H
#define __SRC_PERIODIC_DATA_H

#include <src/math/matrix.h>

namespace bagel {

/* Store data in direct space */
class Data {
  protected:
    int blocksize_;
    int nblock_;

    std::vector<std::shared_ptr<Matrix>> data_;     // (g, i, j)

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & blocksize_ & nblock_ & data_;
    }

  public:
    Data() { }
    Data(const int bsize, const int nblock) : blocksize_(bsize), nblock_(nblock) {
      data_.resize(nblock);
      for (int i = 0; i != nblock; ++i) {
        auto block = std::make_shared<Matrix>(bsize, bsize);
        block->zero();
        data_[i] = block;
      }
    }

    ~Data() { }

    const int blocksize() const { return blocksize_; }
    const int nblock() const { return nblock_; }

    std::shared_ptr<Matrix> operator[] (int i) { assert(i < nblock_ && i >= 0); return data_[i]; };

    std::vector<std::shared_ptr<Matrix>> data() const { return data_; }
    std::shared_ptr<Matrix> data(const int i) const { return data_[i]; }

    void zero() {
      for (auto& block : data_) block->zero();
    }

    void allreduce() {
      for (auto& block : data_) block->allreduce();
    }
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Data)

#endif
