#ifndef __SRC_CI_ZFCI_ZSHCI_H
#define __SRC_CI_ZFCI_ZSHCI_H

#include <src/ci/zfci/zharrison.h>

using namespace std;

namespace bagel {

class ZSHCI : public ZHarrison {
  protected:
    // breit and gaunt
    bool gaunt_;
    bool breit_;
    // To do a one step tight perturbation calculation using Dice.
    bool tight_;
    
    void dump_integrals_and_exit() const override;
    void dump_integrals() const;
    std::shared_ptr<const ZCoeff_Block> init_coeff() override;
  private:
    void write_dice_input();
    void write_guess(ofstream& fs, const int nelea, const int neleb, const int nstate, const int offset);
    void read_dice_rdm();
    void disk_energy();
    void clean();
  private:
    std::shared_ptr<const PTree> shci_data_;
    //specify whether to generate guess by bagel
    bool generate_guess_;
    std::string rdm_file_;
    std::string dice_dir_;
    std::string mpiprefix_;

  public:
    ZSHCI(std::shared_ptr<const PTree> a, std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> b,
           const int ncore = -1, const int nocc = -1, std::shared_ptr<const ZCoeff_Block> coeff_zcas = nullptr, const bool store_c = false, const bool store_g = false);
    void compute() override;
    virtual std::shared_ptr<Kramers<2,ZRDM<1>>> read_external_rdm1(const int ist, const int jst, const std::string& file) const override;
    virtual std::shared_ptr<Kramers<4,ZRDM<2>>> read_external_rdm2(const int ist, const int jst, const std::string& file) const override;
    void compute_rdm12() override;
    void update(std::shared_ptr<const ZCoeff_Block> coeff) override;
    // This function triggers the final step tight calculation of shci.

};

}

#endif
