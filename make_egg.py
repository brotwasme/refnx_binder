import numpy as np
import matplotlib.pyplot as plt
import matplotlib
 
from scipy import special
from refnx.dataset import ReflectDataset
from refnx.analysis import CurveFitter, Objective, Transform, parameter, possibly_create_parameter, Parameters
from refnx.reflect import ReflectModel, SLD, Erf, Component
from refnx.reflect import structure

# |       .        o           .       |
# air-pro, phi-pho, center, mirror (i-o), pro-wat 
#       a-w
#   (x/a)^2 + (y/b)^2 = 1 <-ellipse
# air-pro to cetre = a, b is from area taken up

class param():
    def __init__(real=0, image=0, name='unknown', thickness=0):
        self.real = real
        pass
        

class eggs(Component):#
    def __init__(self, air, hydropho, hydrophil, water, tilt):
        super(eggs)
        names = ['air', 'hydropho', 'hydrophil', 'water']
        more = ['tilt']
        pass

class bsla_thesis(Component): #as in thesis
    def __init__(self, i_as=11.1, sig_as=0.1, i_ap=24.1, th_p= 25.1, i_ps=24.1, n_w=500.1, d_sld=2.01, a=13, b=16.5, c=27.5):
        #interface air solvent, interface width air solvent, interface air protien, interface protien solvent, number of water molecules, sld of deuterated protein 
        super(bsla_thesis).__init__()
        name='bsla'
        self.name = name
        self.interface_air_protein = possibly_create_parameter(0.0,
                name='%s - interface_air_protein' % name)
        self.interface_protein_solvent = possibly_create_parameter(11.1,
                name='%s - interface_protein_solvent' % name)
        self.protein_length = possibly_create_parameter(25.1,
                name='%s - protein_length' % name)
        self.number_of_water_molecules = possibly_create_parameter(500.1,
                name='%s - number_of_water_molecules' % name)
        self.interface_width_air_solvent = possibly_create_parameter(15,
                name='%s - interface_width_air_solvent' % name)
        #self.interface_width_protein_solvent = possibly_create_parameter(0.1,
        #        name='%s - interface_width_protein_solvent' % name)
        self.sld_of_protein = possibly_create_parameter(3.23, #*(10**(-6))
                name='%s - sld_of_protein' % name)
        self.d2o_to_h2o_ratio = possibly_create_parameter(1.0,
                name='%s - d2o_h2o_ratio' % name)
#         self.semi_axis_a = possibly_create_parameter(a,
#                 name='%s - semi_axis_a' % name)
#         self.semi_axis_b = possibly_create_parameter(b,
#                 name='%s - semi_axis_b' % name)
#         self.semi_axis_c = possibly_create_parameter(c,
#                 name='%s - semi_axis_c' % name)
        self.a = 13
        self.b = 16.5
        self.c = 27.5
        self.volume_of_water_molecule = 30
        self.major_axis_length = 55
        self.d = (self.c - self.a)/self.volume_of_water_molecule
        #print([self.interface_air_protein,
        #          self.interface_protein_solvent,
        #          self.protein_length,
        #          self.number_of_water_molecules,
        #          self.interface_width_air_solvent,
        #          #self.interface_width_protein_solvent,
        #          self.sld_of_protein,
        #          self.d2o_to_h2o_ratio])
        
#         b_o=0.5843*10**-4
#         b_d=0.6671*10**-4
#         b_h=-0.3739*10**-4
#         b_d2o = b_o + 2*b_d
#         b_d2o = 1.9185*10**-4
#         b_h2o = -0.1635*10**-4
        
#         self.sld_d2o = b_d2o/self.volume_of_water_molecule

    def __call__(self, z, structure=None):
        self.calculations()
        area_prot = self.area_protein(z)
        area_wat = self.area_water(z)
        #print(area_wat)
        #print('cough')
        #print(area_prot)
        area_total = area_wat+area_prot
        sld_prot = self.sld_protein(area_prot,area_total)
        sld_wat = self.sld_water(area_wat, area_total)
        #print("last",len(sld_prot),len(sld_wat))
        return sld_prot+sld_wat

    @property
    def  parameters(self):
        p = Parameters(name=self.name)
        p.extend([self.interface_air_protein,
                  self.interface_protein_solvent,
                  self.protein_length,
                  self.number_of_water_molecules,
                  self.interface_width_air_solvent,
                  #self.interface_width_protein_solvent,
                  self.sld_of_protein,
                  self.d2o_to_h2o_ratio])
        return p

    def logp(self):
        #full = np.arange(150)
        #a_p = np.array(list(map(self.area_protein, full)))
        #a_w = np.array(list(map(self.area_water, full)))
        #max_w = max(a_w)
        #max_p = max(a_p)
        #max_a = max_w if max_w > max_p else max_p
        return 0

    def calculations(self):
#         a = 13
#         b = 16.5
#         c = 27.5
#         volume_of_water_molecule = 30
        self.e = self.major_axis_length - self.protein_length.value
        self.delta_a = self.a + (self.d*self.e)
        self.delta_c = self.c - (self.d*self.e)

    def area_protein(self, z):
        identity = np.ones(len(z))
        first_bit = np.pi*self.delta_a*self.b/(self.delta_c**2)
        second_bit = z - identity*self.interface_air_protein.value
        final_bit = identity*2*self.delta_c - second_bit
#         print('f',np.pi, self.delta_a,self.b, self.delta_c, first_bit)
#         print(second_bit)
#         print(final_bit)
        #print("pro lengths",len(second_bit), len(final_bit))
        return first_bit*np.multiply(second_bit,final_bit)

    def sld_protein(self, area_p, area_total):
        ratio_area = area_p/area_total
        #print("pro sl lengths",len(area_p),len(ratio_area))
        return ratio_area * self.sld_of_protein.value

    def area_water(self, zs):
        water_length = (self.interface_air_protein.value
                        + self.protein_length.value
                        - self.interface_protein_solvent.value)
        z_s = self.interface_protein_solvent.value+0.5*water_length # calc
        first_bit = (self.number_of_water_molecules.value
                     *self.volume_of_water_molecule/(2*water_length))
        #print(self.interface_width_air_solvent.value, type(self.interface_width_air_solvent.value))
        #print(z,'\n')
        #print(z, z_s ,0.5,self.d, self.interface_width_air_solvent.value,
        #     special.erf((z - z_s + 0.5*self.d)
        #                  /self.interface_width_air_solvent.value))
        second_bits= np.copy(zs)
        #print(z_s, water_length, self.interface_width_air_solvent.value)
        for i in range(len(zs)):
            second_bits[i] = ( special.erf( (zs[i] - z_s + 0.5*water_length)
                              / (self.interface_width_air_solvent.value*np.sqrt(2)) )
                          - special.erf( (zs[i] - z_s - 0.5*water_length)
                                / (self.interface_width_air_solvent.value*np.sqrt(2)) ) )
        #print("wat lengths",len(second_bits))
        return first_bit*second_bits
#             print('helo',
#                   second_bits[i], 
#                   zs[i],
#                   zs[i] - z_s + 0.5*water_length, 
#                   zs[i] - z_s - 0.5*water_length,
#                   special.erf((zs[i] - z_s + 0.5*water_length)
#                     /self.interface_width_air_solvent.value),
#                   special.erf((zs[i] - z_s - 0.5*self.d)
#                     /self.interface_width_air_solvent.value))

    def sld_water(self, area_w, area_total):
        b_d2o = 1.9185#*10**-4
        b_h2o = -0.1635#*10**-4
        sld_d2o = 100*b_d2o/self.volume_of_water_molecule
        sld_h2o = 100*b_h2o/self.volume_of_water_molecule
        sld_solvent = (self.d2o_to_h2o_ratio.value*sld_d2o +
                            (1-self.d2o_to_h2o_ratio.value)*sld_h2o)
        ratio_area = area_w/area_total
        return ratio_area * sld_solvent
    
    def slabs(self,structure=None):
        #in style of https://refnx.readthedocs.io/en/latest/_modules/refnx/reflect/spline.html#Spline
        length = (self.interface_air_protein.value
                        + self.protein_length.value)
        no_slabs = length/1.
        thickness_slabs = length/no_slabs
        slabs = np.zeros((int(no_slabs), 5))
        slabs[:,0] = thickness_slabs
        slabs[-1,3] = 0.5
        distance = np.cumsum(slabs[..., 0]) - 0.5 * thickness_slabs
        slabs[:,1] = self(distance,structure)
        #print("slabs",slabs)
        #slabs  = structure.__profile_slicer(distance,slabs)
        return slabs

#
#air = SLD(value=0+0j, name='air')
#----------
#or 
#re =Parameter(0)
#im=Parameter(0)
#air = SLD([re,im])
#---------
#polymer = SLD(1,'polymer')
#silicon = SLD(2,'silicon')
#structure = air(thick=0,rough=0) | polymer(200,4) | silicon(0,3)
# air-polymer roughness of 4, polymer size of 200
#  #Erf()  <-error function
# structure[1].interfaces = Erf() # air-polymer interface
# structure[2].interfaces = Erf()
#need lnprob
#sort vol fractions 
#add areas etc.
#possibly using possibly` creat parameter and parameter p.extend(
#tilt and other values not used either


class bsla_non(Component): #not as in thesis
    def __init__(self, area_per_molecule=1, angle=1, length_of_molecule=1):
        #area per molecule, 
        super(bsla).__init__()
        name='bsla'
        pass

    
    #for calculating unphysical thing
    def volume_fraction_protien(self, z):
        pass
    
    def volume_fraction_air(self, z):
        pass
    
    def volume_fraction_water(self, z):
        pass
    
    
class bsla_equation(Component): # equation based
    def __init__(self, extent, vs, dz,name='bsla', microslab_max_thickness=1):
        pass
    
    def __call__(self, z):
        pass
    
def test_run():
    print("\n\n\n")
    from refnx.dataset import Data1D
    from refnx.dataset import ReflectDataset
    import refnx
    import data_in
    data = data_in.data_in('d2o/29553_54.dat')
    # dataset = data # ...
    data = Data1D(data) 
    from make_egg import bsla_thesis
#     air = SLD(0)
#     air = air(0,0)
    bt = bsla_thesis()
    bt.interface_protein_solvent.setp(vary=True, bounds=(11, 40))
    bt.protein_length.setp(vary=True, bounds=(25, 55))
    bt.number_of_water_molecules.setp(vary=True, bounds=(1, 10000))
    bt.interface_width_air_solvent.setp(vary=True, bounds=(0.001, 30))
    #bt.interface_width_protein_solvent.setp(vary=True, bounds=(0, 5))
    bt.sld_of_protein.setp(vary=True, bounds=(1.92, 6.21)) # *(10**(-6))
    bt.d2o_to_h2o_ratio.setp(vary=True, bounds=(0, 1))
#     if isinstance(bt, Component):
#         print("it is comp")
#     if isinstance(bt, Structure):
#         print("it is")
    #print(bt.parameters)
#     d2o = 1.9185/0.3
#     h2o = -0.1635/0.3
#     solvent = SLD((bt.d2o_to_h2o_ratio.value*d2o + (1-bt.d2o_to_h2o_ratio.value)*h2o))
#     solvent = solvent(0,0)
#     from refnx.reflect import Structure
#     structure = air|bt|solvent
#     structure.name = "bsla"
    from refnx.reflect import ReflectModel
    model = ReflectModel(bt)#structure)
    from refnx.analysis import Transform, CurveFitter, Objective
    objective = Objective(model,data)
    fitter = CurveFitter(objective)
    fitter.fit('differential_evolution');
    import matplotlib.pyplot as plt
    #%matplotlib notebook
#     plt.plot(*bt.sld_profile())
    objective.plot()
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Q')
    plt.ylabel('Reflectivity')
    plt.legend()
    print(bt)
    plt.show()
    
if __name__ == '__main__':
    test_run()