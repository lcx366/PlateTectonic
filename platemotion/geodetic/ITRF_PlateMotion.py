from .sites_select import sites_plate
from .estimate_eulerpole import PMC_iterate

def ITRF_PlateMotion(modelname,dir_to=None,out_days=90):

    platemodel = sites_plate(modelname,dir_to,out_days)
    platemodel_info = vars(platemodel)
    print('Estimating Euler Pole ... ',end='')
    for platename in platemodel_info.keys():
        plate = platemodel_info[platename]
        sites_info = plate.sites
        epr,sites_retain_info = PMC_iterate(sites_info)
        plate.set_epr(epr)
        plate.set_omega(epr['omega_cartesian'],'cartesian')
        plate.set_sites(sites_retain_info)
        platemodel.add_plate(plate)
    print('Finished')    
    return platemodel