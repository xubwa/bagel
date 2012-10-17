//
// BAGEL - Parallel electron correlation program.
// Filename: test_molden.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker < shane.parker@u.northwestern.edu >
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#include <src/molden/molden.h>
#include <src/scf/fock.h>

double molden_out_energy(std::string inp1, std::string inp2) {

    std::shared_ptr<std::ofstream> ofs(new std::ofstream(inp1 + ".testout", std::ios::trunc));
    std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

    // a bit ugly to hardwire an input file, but anyway...
  {
    std::shared_ptr<InputData> idata(new InputData("../../test/" + inp1 + ".in"));
    std::shared_ptr<Geometry> geom(new Geometry(idata->get_input("molecule")));
    std::list<std::pair<std::string, std::multimap<std::string, std::string> > > keys = idata->data();

    std::shared_ptr<Reference> ref;

    for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
      if (iter->first == "df-hf") {
        std::shared_ptr<SCF<1> > scf(new SCF<1>(iter->second, geom));
        scf->compute();
        ref = scf->conv_to_ref();
      }
      else if (iter->first == "print") {
        std::multimap<std::string, std::string> pdata = iter->second;
        bool orbitals = read_input<bool>(pdata, "orbitals", false);
        std::string out_file = read_input<std::string>(pdata, "file", "test.molden");
     
        Molden molden(geom->spherical());
        molden.write_geo(geom, out_file);
        if (orbitals) molden.write_mos(ref, out_file);

      }
    }
  }

  // a bit ugly to hardwire an input file, but anyway...
  {
    std::shared_ptr<InputData> idata(new InputData("../../test/" + inp2 + ".in"));
    std::shared_ptr<Geometry> geom(new Geometry(idata->get_input("molecule")));
    std::list<std::pair<std::string, std::multimap<std::string, std::string> > > keys = idata->data();

    Molden mfs;

    std::shared_ptr<const Coeff> coeff = mfs.read_mos(geom, "test.molden");

    std::shared_ptr<Matrix1e> ao_density = coeff->form_density_rhf(geom->nele()/2);
    std::shared_ptr<Fock<1> > hcore(new Fock<1>(geom));
    std::shared_ptr<Fock<1> > fock(new Fock<1>(geom, hcore, ao_density, geom->schwarz()));

    Matrix1e hcore_fock = (*hcore + *fock);
    double energy = ((*ao_density)*(hcore_fock)).trace();//->ddot(*hcore_fock.transpose());
    energy = 0.5*energy + geom->nuclear_repulsion();

    std::cout << energy << std::endl;
    std::cout << geom->nele() << std::endl;

    return energy;
  }
}

BOOST_AUTO_TEST_SUITE(TEST_MOLDEN)

BOOST_AUTO_TEST_CASE(MOLDEN_OUTIN) {
    BOOST_CHECK(compare(molden_out_energy("hf_writemolden", "hf_readmolden"),-99.84772354 ));
}

BOOST_AUTO_TEST_SUITE_END()
