
import math
from collections import OrderedDict
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import \
    Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from gen_helper import isLastCopy
from utils import deltaR


class GenParticleModule(Module):

    def __init__(
        self,
        inputGenCollection=lambda event: Collection(event, "GenPart"),
        outputName="genPart",
        storeKinematics= ['pt','eta','phi','mass'],
    ):
        self.inputGenCollection = inputGenCollection
        self.outputName = outputName
        self.storeKinematics = storeKinematics

        self.verbose = False

        self.quarks_list = [1,2,3,4]
        self.leptons_list = [11,13,15]

        self.genTopKeys = [
            'first_daughter', 'second_daughter', 'third_daughter', 'has_hadronically_decay',]
        if Module.globalOptions['isSignal']: self.genTopKeys.append('from_resonance')

    def has_x_as_mother(self, genParticle_idx, target_mother_pdg, event):
        # recursive check if x is the mother of the particles under investigation
        gen_collection = self.inputGenCollection(event)

        while genParticle_idx != -1:
            mother_idx = gen_collection[genParticle_idx].genPartIdxMother

            if mother_idx != -1:
                mother_pdg = gen_collection[mother_idx].pdgId
                if abs(mother_pdg) == target_mother_pdg:
                    return True
                genParticle_idx = mother_idx
            else:
                break
        return False

    def has_x_as_daughter(self, genP, daughter_of_interest, event):
        # recursively check if daughter with pdgId x is a daughter
        genparticles = map(lambda genP: genP, self.inputGenCollection(event))[genP._index+1:]

        if not genparticles:
            if self.verbose:
                print('No daughter of pdg id {} found for this genP'.format(daughter_of_interest))
            return False, None

        genp_mothers = [genp.genPartIdxMother for genp in genparticles]

        for igenP, genparticle in enumerate(genparticles):
            is_last_particle = igenP + 1 == len(genparticles) # flag if the genp under investigation is the last one in the list
            is_daughter = genP._index == genparticle.genPartIdxMother and genparticle.pdgId == daughter_of_interest

            if is_daughter:
                if self.verbose:
                    print('The {} in pos {} has a daughter {} in pos {}'.format(genP.pdgId, genP._index, genparticle.pdgId, genparticle._index))

                is_last_copy = isLastCopy(genparticle) or genparticle._index not in genp_mothers
                # additional check with "genparticle._index not in genp_mothers":
                # it can happens that a genparticle is not considered last copy even
                # if it has no daughter (e.g. eventual daughters are not recorded due to genfilters)

                if is_last_copy or is_last_particle:
                    if self.verbose: print('which is the last copy')
                    return True, genparticle
                else:
                    if self.verbose: print('which is not the last copy.')
                    genP = genparticle
                    return self.has_x_as_daughter(genP, daughter_of_interest, event)

        return False, None

    def check_W_daughters(self, W, daughters_list, event):
        W_sign = np.sign(W.pdgId)

        if W_sign > 0: 
            decay_pairs = [
                (2, -1), (2, -3), (2, -5),
                (4, -1), (4, -3), (4, -5),
                (6, -1), (6, -3), (6, -5),
                (-11, 12), (-13, 14), (-15, 16)
            ]
        else:
            decay_pairs = [
                (-2, 1), (-2, 3), (-2, 5),
                (-4, 1), (-4, 3), (-4, 5),
                (-6, 1), (-6, 3), (-6, 5),
                (11, -12), (13, -14), (15, -16)
            ]

        # recursively check the W daughters
        for daughter_1_id, daughter_2_id in decay_pairs:
            check_1, daughter_1 = self.has_x_as_daughter(W, daughter_1_id, event)
            check_2, daughter_2 = self.has_x_as_daughter(W, daughter_2_id, event)
            if check_1 and check_2:
                daughters_list.extend([daughter_1, daughter_2])
                break

    def deltaPhi(self,phi1, phi2):
        res = phi1-phi2
        while (res > math.pi):
            res -= 2 * math.pi
        while (res <= -math.pi):
            res += 2 * math.pi

        return res

    def relDeltaPt(self,pt1, pt2):
        if pt2 != 0.0:
            return abs((pt1 - pt2) / pt2)
        return 1.

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("ngenTops", "I")
        for genTopKey in self.genTopKeys:
            if 'deltaR' in genTopKey or 'rho_over_pt' in genTopKey: 
                self.out.branch("genTops_"+genTopKey, "F", lenVar="ngenTops")
            else: 
                self.out.branch("genTops_"+genTopKey, "I", lenVar="ngenTops")
        for variable in self.storeKinematics:
            self.out.branch("genTops_"+variable, "F", lenVar="ngenTops")
        # self.out.branch("ngenTop_daughters","I")
        self.out.branch("genTops_daughters_pdgId", "I", lenVar="ngenTops")

        for genlepton in ['genElectrons', 'genMuons']:
            self.out.branch("n{}".format(genlepton), "I")
            for variable in self.storeKinematics:
                self.out.branch("{}_{}".format(genlepton, variable), "F", lenVar="n{}".format(genlepton))
            self.out.branch("{}_status".format(genlepton), "F", lenVar="n{}".format(genlepton))
            self.out.branch("{}_index".format(genlepton), "F", lenVar="n{}".format(genlepton))

        # self.out.branch("ngenParticle", "I")
        # for variable in self.storeKinematics:
        #     self.out.branch("genParticle_"+variable, "F", lenVar="ngenParticle")  
        # for genvar in ['pdgId', 'genPartIdxMother', 'status', 'statusFlags']:
        #     self.out.branch("genParticle_"+genvar, "I", lenVar="ngenParticle")  


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        genParticles = self.inputGenCollection(event)

        if self.verbose:
            print('########## EVENT ############')
            pdgIds = [getattr(genp, "pdgId") for genp in genParticles]
            genPartIdxMothers = [getattr(genp, "genPartIdxMother") for genp in genParticles]

            print("## pdgIds")
            print(pdgIds)
            print("## genPartIdxMother")
            print(genPartIdxMothers)

            print("## positional pdgIds")
            print({num: genp.pdgId for num, genp in enumerate(genParticles)})
            print("## positional genPartIdxMother")
            print({num: genp.genPartIdxMother for num, genp in enumerate(genParticles)})

        # ---- creating new list of last copies of t, W, b, resonance
        # these list are passed to another module to calculated the decays products of t, W, b which are inside jets (for jet-classification)
        quarks_from_w = [1, 2, 3, 4]
        gen_top_quarks = []
        gen_w_bosons_not_from_top_list, gen_b_quarks_not_from_top = [], []
        gen_particles_not_from_top = []
        gen_electrons = []
        gen_muons = []

        for genParticle in genParticles:
            if abs(genParticle.pdgId) == 11:# & genParticle.status == 1: 
                gen_electrons.append(genParticle)
            if abs(genParticle.pdgId) == 13:# & genParticle.status == 1: 
                gen_muons.append(genParticle)

            if isLastCopy(genParticle):
                if abs(genParticle.pdgId) == 6:
                    if self.verbose:
                        print('GenTop last copy in pos {}'.format(genParticle._index))

                    # ---- check the origin of the genTop (if comes from resonance)
                    if Module.globalOptions['isSignal']:
                        setattr(genParticle, 'from_resonance', False)
                        if self.has_x_as_mother(genParticle._index, 6000046, event): #6000055 (Z') 6000045 (tth) 6000046 (tta) 
                            setattr(genParticle, 'from_resonance', True)
                            if self.verbose:
                                print('...which comes from resonance!')
                                print(genParticle.from_resonance)
                    # ----

                    # ---- genTop daughter search
                    # the function "has_x_as_daughter" recursively finds the last copy
                    # of the daughter x whose pdgId is given as second argument
                    gentop_daughters = []
                    b_as_daughter = self.has_x_as_daughter(genParticle, np.sign(genParticle.pdgId)*5, event)
                    w_as_daughter = self.has_x_as_daughter(genParticle, np.sign(genParticle.pdgId)*24, event)

                    if b_as_daughter[0]:
                        gentop_daughters.append(b_as_daughter[1])

                    if w_as_daughter[0]:
                        self.check_W_daughters(w_as_daughter[1], gentop_daughters, event)

                    setattr(genParticle, 'daughters', [])
                    for n_daugther in ['first_daughter','second_daughter','third_daughter']:
                        setattr(genParticle, n_daugther, -99) # attribute initialization
                    setattr(genParticle, 'daughters', gentop_daughters)
                    if len(gentop_daughters)==len(['first_daughter','second_daughter','third_daughter']):
                        for i, daughter in zip(['first_daughter','second_daughter','third_daughter'], gentop_daughters):
                            setattr(genParticle, i, daughter.pdgId)
                    else:
                        if self.verbose:
                            print(gentop_daughters)
                            print("Wierd Event...")
                    if self.verbose:
                        print('Top pos.[{}] has {} in pos.[{}] as daughters'.format(genParticle._index, list(map(lambda daughter: daughter.pdgId, gentop_daughters)), list(map(lambda daughter: daughter._index, gentop_daughters))))

                    setattr(genParticle, 'has_hadronically_decay', False)
                    if any(top_daughter in list(map(lambda daughter: daughter.pdgId, gentop_daughters)) for top_daughter in quarks_from_w):
                        if self.verbose: print('Top pos.[{}] decays hadronically')
                        setattr(genParticle, 'has_hadronically_decay', True)
                    # ----
                    gen_top_quarks.append(genParticle)

                elif abs(genParticle.pdgId)==5 and not self.has_x_as_mother(genParticle._index, 6, event):
                    gen_b_quarks_not_from_top.append(genParticle)

                elif abs(genParticle.pdgId)==24 and not self.has_x_as_mother(genParticle._index, 6, event):
                    gen_w_bosons_not_from_top_list.append(genParticle)

                    genW_daughters = [] #{'from_resonance': [], 'not_from_resonance': []}
                    self.check_W_daughters(genParticle, genW_daughters, event)

                    setattr(genParticle, 'daughters', [])
                    if len(genW_daughters)==len(['first_daughter','second_daughter']):
                        setattr(genParticle, 'daughters', genW_daughters)
                    else:
                        if self.verbose:
                            print("Wierd Event...")

                else:
                    mother_24 = self.has_x_as_mother(genParticle._index, 24, event)
                    mother_6 = self.has_x_as_mother(genParticle._index, 6, event)
                    mother_5 = self.has_x_as_mother(genParticle._index, 5, event)

                    leptons_set = set(self.leptons_list)
                    quarks_set = set(self.quarks_list)

                    if (not (mother_24 or mother_6) or not (mother_5 or mother_6)) and (abs(genParticle.pdgId) in leptons_set or abs(genParticle.pdgId) in quarks_set):
                        gen_particles_not_from_top.append(genParticle)

            else: continue
        # ----

        setattr(event, 'genTops', gen_top_quarks)
        setattr(event, 'gen_w_bosons_not_from_top', gen_w_bosons_not_from_top_list)
        setattr(event, 'gen_b_quarks_not_from_top', gen_b_quarks_not_from_top)
        setattr(event, 'gen_particles_not_from_top', gen_particles_not_from_top)

        # self.out.fillBranch("ngenParticle", len(genParticles))
        # for variable in self.storeKinematics:
        #     self.out.fillBranch("genParticle_"+variable, map(lambda genp: getattr(genp, variable), genParticles))  
        # for genvar in ['pdgId', 'genPartIdxMother', 'status', 'statusFlags']:
        #     self.out.fillBranch("genParticle_"+genvar, map(lambda genp: getattr(genp, genvar), genParticles))  

        self.out.fillBranch("ngenTops", len(gen_top_quarks))
        for genTopKey in self.genTopKeys:
            self.out.fillBranch("genTops_"+genTopKey, map(lambda gentop: getattr(gentop,genTopKey), gen_top_quarks))
        for variable in self.storeKinematics:
            self.out.fillBranch("genTops_"+variable, map(lambda gentop: getattr(gentop,variable), gen_top_quarks))

        for genlepton, genlepton_collection in zip(['genElectrons', 'genMuons'], [gen_electrons, gen_muons]):
            self.out.fillBranch("n{}".format(genlepton), len(genlepton_collection))
            for variable in self.storeKinematics:
                self.out.fillBranch("{}_{}".format(genlepton, variable), map(lambda genlep: getattr(genlep, variable), genlepton_collection))
            self.out.fillBranch("{}_status".format(genlepton), map(lambda genlep: getattr(genlep, 'status'), genlepton_collection))
            self.out.fillBranch("{}_index".format(genlepton), map(lambda genlep: getattr(genlep, '_index'), genlepton_collection))

        return True

