import sys, unittest
from md import calcenergy
from ase.lattice.cubic import FaceCenteredCubic
from asap3 import EMT

class MdTests(unittest.TestCase):
      
      def test_calcenergy(self):
            # Set up a crystal
            atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                              symbol="Cu",
                              size=(10, 10, 10),
                              pbc=True)
            atoms.calc = EMT()
            res = calcenergy(atoms)

            if len(res) == 4:
                  self.assertTrue(True)
            else:
                  self.assertTrue(False)
         
          
if __name__ == "__main__":
   tests = [unittest.TestLoader().loadTestsFromTestCase(MdTests)]
   testsuite = unittest.TestSuite(tests)
   result = unittest.TextTestRunner(verbosity=0).run(testsuite)
   sys.exit(not result.wasSuccessful())


