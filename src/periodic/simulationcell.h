//
// BAGEL - Parallel electron correlation program.
// Filename: simulationcell.h
// Copyright (C) 2015 Toru Shiozaki
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


#ifndef __BAGEL_SRC_PERIODIC_SIMULATIONCELL_H
#define __BAGEL_SRC_PERIODIC_SIMULATIONCELL_H

#include <src/wfn/geometry.h>

namespace bagel {

class SimulationCell { /* cubic, same or larger than primitive cell */
  protected:
    std::vector<std::shared_ptr<const Atom>> atoms_;
    int ndim_;
    std::vector<std::array<double, 3>> primitive_vectors_; // orthogonal vectors
    std::array<double, 3> charge_centre_;

    int num_jvectors_;
    std::vector<std::array<double, 3>> jvectors_;
    int nbasis_;
    void init();
    int ws_;
    double extent_, radius_;
    std::vector<std::shared_ptr<const ZMatrix>> multipoles_;
    void compute_multipoles(const int lmax);
    void compute_extent(const double thresh = PRIM_SCREEN_THRESH);

  public:
    SimulationCell() { }
    SimulationCell(const std::shared_ptr<const Geometry> geom);
    SimulationCell(std::vector<std::shared_ptr<const Atom>> atoms, const int ndim, std::vector<std::array<double, 3>> prim_vec);
    virtual ~SimulationCell() { }

    int ndim() const { return ndim_; }
    int num_jvectors() const { return num_jvectors_; }

    std::vector<std::array<double, 3>> jvectors() const { return jvectors_; }
    std::array<double, 3> jvectors(const int i) const { return jvectors_[i]; }

    std::vector<std::shared_ptr<const ZMatrix>> multipoles() const { return multipoles_; }

    std::array<double, 3> centre() const { return charge_centre_; }
    double centre(const int i) const { return charge_centre_[i]; }
    double extent() const { return extent_; }
    double radius() const { return radius_; }
};

}

#endif
