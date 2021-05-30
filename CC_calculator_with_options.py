import numpy as np
from free_energy_codes.component_contribution.component_contribution import *
from free_energy_codes.component_contribution.thermodynamic_constants import default_RT
# from free_energy_codes.component_contribution.component_contribution_trainer import ComponentContribution
# from component_contribution.kegg_model import KeggModel
from free_energy_codes.component_contribution.compound import Compound
import os
import pdb
import pandas as pd

class dGReaction:
    def __init__(self, S, cids, inchis,rids=None):
        self.S = S
        self.cids = cids
        self.inchis = inchis
        assert len(self.cids) == self.S.shape[0]

        # remove H+ from the stoichiometric matrix if it exists
        if 'C00080' in self.cids:
            i = self.cids.index('C00080')
            self.S = np.vstack((self.S[:i,:], self.S[i+1:,:]))
            self.cids.pop(i)

    def get_dG0(self,cc,options):
        """
            returns the estimated dG0
        """
        dG0_compounds = np.matrix(np.zeros((self.S.shape[0], 1)))
        for i, cid in enumerate(self.cids):
            if options[i] == 1:
                dG0_compounds[i, 0] = cc.get_major_ms_dG0_f(cid)
            if options[i] == 2:
                if "InChI=" not in self.inchis[i]:
                    print 'invalid inchi string'
                    break
                dG0_compounds[i, 0] = cc.get_major_ms_dG0_f_from_InChI(self.inchis[i], database='Test', compoundID=cid)
        print dG0_compounds

        dG0 = np.dot(self.S.T, dG0_compounds)
        return dG0

    def get_transformed_dG0(self, dG0, pH, I, T,options):
        """
            returns the estimated dG0_prime
        """
        dG0_prime = dG0 + self.get_transform_ddG0(pH=pH, I=I, T=T,options=options)
        return dG0_prime


    def get_transform_ddG0(self,pH, I, T,options):
        """
        needed in order to calculate the transformed Gibbs energies of the
        model reactions.

        Returns:
            an array (whose length is self.S.shape[1]) with the differences
            between DrG0_prime and DrG0. Therefore, one must add this array
            to the chemical Gibbs energies of reaction (DrG0) to get the
            transformed values
        """
        ddG0_compounds = np.matrix(np.zeros((self.S.shape[0], 1)))
        for i, cid in enumerate(self.cids):
            # comp = self.ccache.get_compound(cid)
            if options[i] == 1:
                comp = Compound.from_kegg(cid)
            if options[i] == 2:
                if "InChI=" not in self.inchis[i]:
                    print 'invalid inchi string'
                    break
                comp = Compound.from_inchi('Test', cid, self.inchis[i])

            ddG0_compounds[i, 0] = comp.transform_pH7(pH, I, T)

        ddG0_forward = np.dot(self.S.T, ddG0_compounds)
        return ddG0_forward

    def get_dGm(self,dG0_prime):
        """
        # reclculate dG with concentrations= 1mM
        """
        mM_conc = 1e-3 * np.matrix(np.ones((len(self.cids), 1)))
        if 'C00001' in self.cids:
            mM_conc[self.cids.index('C00001'), 0] = 1.0 # concentarion of water is ignored
        dGm_prime = dG0_prime + default_RT * self.S.T * np.log(mM_conc)
        return dGm_prime

if __name__ == '__main__':
    ### Example: calculate dG, dG'0, and dG'm for reaction R00258
    # - dG is the standard Gibss free energy without considering pH, ionic strength
    # - dG'0 is the change in Gibbs free energy due to a chemical reaction due a reaction at a particular pH and ionic strength
    # - dG'm is for much more appropriate standard concentration (1 mM)
    # - R00258: L-Alanine + 2-Oxoglutarate <=> Pyruvate + L-Glutamate
    #           C00041    + C00026         <=> C00022   + C00025

    # input information
    df = pd.read_csv('dG_with_options_input.csv',index_col=0)
    S = np.asarray(df['stoich'].tolist())
    cids = df['cid'].tolist()
    inchis = df['inchi'].tolist()
    options = df['options'].tolist()  # how to calculate dG, 1 = use cids (kegg id) 2= use inChi (for novel metabolites)
    
    # start component contribution
    cc = ComponentContribution.init()
    reaction = dGReaction(S,cids,inchis)

    # calculate gibbs free energy for different options
    dG0 = reaction.get_dG0(cc,options) # standard
    dG0_prime = reaction.get_transformed_dG0(dG0,7.0, 0.1, 298.15,options) # pH = 7.0, I = 0.1, and T = 298.15
    dGm_prime = reaction.get_dGm(dG0_prime) # concentration = 1mM

    print "dG0  = %8.1f" % (dG0) + ' kJ/mol'
    print "dG'0 = %8.1f" % (dG0_prime) + ' kJ/mol'
    print "dG'm = %8.1f" % (dGm_prime) + ' kJ/mol'