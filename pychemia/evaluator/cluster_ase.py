import os
import time
from multiprocessing import Pool, Process
import pychemia
from pychemia.utils.serializer import generic_serializer
from pychemia.external.ase import pychemia2ase,ase2pychemia
from pychemia.external.pymatgen import pychemia2pymatgen
from ase.optimize import BFGS,QuasiNewton
from ase.io import read
from deepmd.calculator import DP
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import PointGroupAnalyzer

__author__ = 'Haidi Wang'

cal=DP(model="frozen_model_cu.pb",type_dict={'Cu':0})
aaa=AseAtomsAdaptor()
print(cal)

def cluster_ase_worker(db_settings):

    while True:
        pcdb = pychemia.db.get_database(db_settings)
        population = pychemia.population.LJCluster(pcdb)

        entry = population.pcdb.db.pychemia_entries.find_one({'status.' + population.tag: True,
                                                              'status.lock': {'$exists': False},
                                                              'properties': {}}, {'_id': 1})

        print(entry)
        if entry is not None:
            population.pcdb.lock(entry['_id'])
            structure = population.pcdb.get_structure(entry['_id'])
            atoms=pychemia2ase(structure)
            max_distance=atoms.get_all_distances().max()
            atoms.set_cell([max_distance+15]*3)
            atoms.set_pbc([True,True,True])
            atoms.center()
            atoms.set_calculator(cal)
            dyn = BFGS(atoms,logfile='tmp.log')
            ret=dyn.run(fmax=1e-3,steps=5000)
            if ret:
               print("converge in %d steps"%dyn.nsteps)
               energy=atoms.get_potential_energy()[0]
            else:
               print("unconverge")
               energy='622427' 
            forces=generic_serializer(atoms.get_forces())
            atoms.set_pbc(False)
            molecule=aaa.get_molecule(atoms)
            pga=PointGroupAnalyzer(molecule)
            pg=str(pga.get_pointgroup())
            print("----->%s"%pg)
            structure=ase2pychemia(atoms)
            properties = {'forces': forces, 'energy': energy ,'point_group':pg}
            population.pcdb.update(entry['_id'], structure=structure, properties=properties)
            population.pcdb.unlock(entry['_id'])
        else:
            break


def cluster_ase_evaluator(db_settings, nparal):
    pcdb = pychemia.db.get_database(db_settings)
    population = pychemia.population.LJCluster(pcdb)
    print('Staring evaluator for ', population.name)
    while True:
        entry = population.pcdb.db.pychemia_entries.find_one({'status.' + population.tag: True,
                                                              'status.lock': {'$exists': False},
                                                              'properties': {}}, {'_id': 1})

        if entry is None:
            time.sleep(2)
            create_pool = False
        else:
            create_pool = True

        print("create_pool %s"%create_pool)
        if create_pool:
            pool = Pool(processes=nparal)
            pool.map(cluster_ase_worker, nparal * [db_settings])
            pool.close()
            pool.join()


def cluster_ase_launcher(db_settings, nparal):
    p = Process(target=cluster_ase_evaluator, args=(db_settings, nparal))
    p.start()
    return p

def main(db_settings):
    pcdb = pychemia.db.get_database(db_settings)
    population = pychemia.population.LJCluster(pcdb)
    print('Staring evaluator for ', population.name)
    while True:
        entry = population.pcdb.db.pychemia_entries.find_one({'status.' + population.tag: True,
                                                              'status.lock': {'$exists': False},
                                                              'properties': {}}, {'_id': 1})

        if entry is None:
            time.sleep(2)
            create_pool = False
        else:
            create_pool = True

        print("create_pool %s"%create_pool)
        cluster_ase_worker(db_settings)

if __name__=='__main__':
   from monty.serialization import dumpfn,loadfn
   main(loadfn('db_settings.json'))
