import unittest
from tRNA_position import *

class TestSSAlignment(unittest.TestCase):
  def setUp(self):
    # hg19-tRNA-His-GTG-1-2, aligned on numeric model
    self.positions = annotate_positions('(((..(((.(......,,<<<<_............_______..............___>>>>,<<<<<______......................................._>.>>>>,.,<<<<<<<_____>>>>>>.>,,,<<<<<_______>>>>.>.......)))))))::::')

  def test_anticodon_with_gaps_on_5prime_side(self):
    correct_positions = zip([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], ['1:179', '2:178', '3:177', '4', '5', '6:176', '7:175', '8:174', '9', '10:173'])
    correct_sprinzl = zip([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], ['1:72', '2:71', '3:70', '3i1', '3i2', '4:69', '5:68', '6:67', '6i1', '7:66'])
    for index, position in correct_positions:
      self.assertEqual(self.positions[index].position, position, 'wrong position, got {} but expected {}'.format(self.positions[index].position, position))
    for index, sprinzl in correct_sprinzl:
      self.assertEqual(self.positions[index].sprinzl, sprinzl, 'wrong position, got {} but expected {}'.format(self.positions[index].sprinzl, sprinzl))

  def test_dstem_gapless(self):
    correct_positions = zip([18, 19, 20, 21], ['19:63', '20:62', '21:61', '22:60'])
    correct_sprinzl = zip([18, 19, 20, 21], ['10:25', '11:24', '12:23', '13:22'])
    for index, position in correct_positions:
      self.assertEqual(self.positions[index].position, position, 'wrong position, got {} but expected {}'.format(self.positions[index].position, position))
    for index, sprinzl in correct_sprinzl:
      self.assertEqual(self.positions[index].sprinzl, sprinzl, 'wrong position, got {} but expected {}'.format(self.positions[index].sprinzl, sprinzl))

  def test_tstem_with_gap_on_3prime_side(self):
    correct_positions = zip([131, 132, 133, 134, 135], ['148:165', '149:163', '150:162', '151:161', '152:160'])
    correct_sprinzl = zip([131, 132, 133, 134, 135], ['49:65', '50:64', '51:63', '52:62', '53:61'])
    for index, position in correct_positions:
      self.assertEqual(self.positions[index].position, position, 'wrong position, got {} but expected {}'.format(self.positions[index].position, position))
    for index, sprinzl in correct_sprinzl:
      self.assertEqual(self.positions[index].sprinzl, sprinzl, 'wrong position, got {} but expected {}'.format(self.positions[index].sprinzl, sprinzl))


if __name__ == '__main__':
    unittest.main()