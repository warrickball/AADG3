import numpy as np
import AADG3
import unittest

tmpfile = 'tmpfile~'

class TestAADG3(unittest.TestCase):
    
    def test_namelist_io(self):
        nml0 = AADG3.load_namelist('test.in')
        self.assertEqual(nml0['user_seed'], 1000)
        self.assertEqual(nml0['n_relax'], 2000)
        self.assertEqual(nml0['n_cadences'], 3000)
        self.assertEqual(nml0['n_fine'], 50)
        self.assertEqual(nml0['cadence'], 60e0)
        self.assertEqual(nml0['sig'], 250e0)
        self.assertEqual(nml0['rho'], 0e0)
        self.assertEqual(nml0['tau'], 120e0)
        self.assertEqual(nml0['p(1)'], 1e0)
        self.assertEqual(nml0['p(2)'], 1e0)
        self.assertEqual(nml0['p(3)'], 1e0)
        self.assertEqual(nml0['inclination'], 45e0)
        self.assertEqual(nml0['cycle_period'], 100e0)
        self.assertEqual(nml0['cycle_phase'], 0e0)
        self.assertEqual(nml0['add_granulation'], True)
        self.assertEqual(nml0['modes_filename'], 'basic.con')
        self.assertEqual(nml0['rotation_filename'], 'basic.rot')
        self.assertEqual(nml0['output_filename'], 'basic.out')

        AADG3.save_namelist(tmpfile, nml0)
        nml1 = AADG3.load_namelist(tmpfile)

        for k in nml0:
            self.assertEqual(nml0[k], nml1[k])

        for k in nml1:
            self.assertEqual(nml0[k], nml1[k])

    def test_modes_io(self):
        modes0 = AADG3.load_modes('test.con')

        # ('l', 'n', 'freq', 'width', 'power', 'dfreq')
        for row in modes0:
            self.assertEqual(row['n'], 21)
            self.assertEqual(row['dfreq'], 0.0)

        np.savetxt(tmpfile, modes0, fmt='%4i%4i%16.7f%16.7f%16.7f%16.7f')
        modes1 = AADG3.load_modes(tmpfile)

        self.assertTrue(np.all(modes0 == modes1))


    def test_rot_io(self):
        rot0 = AADG3.load_rot('test.rot')

        for row in rot0:
            self.assertEqual(row['n'], 21)
            self.assertTrue(row['m'] <= row['l'])
            self.assertTrue(row['m'] > 0)

        np.savetxt(tmpfile, rot0, fmt='%4i%4i%4i%16.7f')
        rot1 = AADG3.load_rot(tmpfile)

        self.assertTrue(np.all(rot0 == rot1))


    def test_generate_const_rot(self):
        splitting = np.random.rand()
        modes = AADG3.load_modes('test.con')

        rot = AADG3.generate_const_rot(modes, splitting=splitting)

        for row in rot:
            self.assertEqual(row['splitting'], splitting)
            self.assertEqual(row['n'], 21)
            self.assertTrue(row['m'] <= row['l'])
            self.assertTrue(row['m'] > 0)
        
        
if __name__ == '__main__':
    unittest.main()
