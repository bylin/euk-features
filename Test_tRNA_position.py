import unittest
from nose.tools import nottest
from tRNA_position import *

# Run with: nosetests Test_tRNA_position.py

class AlignmentTestCase(unittest.TestCase):
  @nottest
  def testPositions(self, indices, positions, sprinzl):
    for index, position in zip(indices, positions):
      self.assertEqual(self.positions[index].position, position, 'wrong position, got {} but expected {}'.format(self.positions[index].position, position))
    for index, sprinzl in zip(indices, sprinzl):
      self.assertEqual(self.positions[index].sprinzl, sprinzl, 'wrong position, got {} but expected {}'.format(self.positions[index].sprinzl, sprinzl))

class TestSSAlignmentHis(AlignmentTestCase):
  @classmethod
  def setUpClass(self):
    # hg19-tRNA-His-GTG-1-2, aligned on numeric model
    self.positions = annotate_positions('(((..(((.(......,,<<<<_............_______..............___>>>>,<<<<<______......................................._>.>>>>,.,<<<<<<<_____>>>>>>.>,,,<<<<<_______>>>>.>.......)))))))::::')

  def test_acceptor_with_gaps_on_5prime_side_plus_inserts(self):
    indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    positions = ['1:179', '2:178', '3:177', '4', '5', '6:176', '7:175', '8:174', '9', '10:173']
    sprinzl = ['1:72', '2:71', '3:70', '3i1', '3i2', '4:69', '5:68', '6:67', '6i1', '7:66']
    self.testPositions(indices, positions, sprinzl)

  def test_dstem_gapless(self):
    indices = [18, 19, 20, 21]
    positions = ['19:63', '20:62', '21:61', '22:60']
    sprinzl = ['10:25', '11:24', '12:23', '13:22']
    self.testPositions(indices, positions, sprinzl)

  def test_tstem_with_gap_on_3prime_side(self):
    indices = [131, 132, 133, 134, 135]
    positions = ['148:165', '149:163', '150:162', '151:161', '152:160']
    sprinzl = ['49:65', '50:64', '51:63', '52:62', '53:61']
    self.testPositions(indices, positions, sprinzl)

class TestSSAlignmentAla(AlignmentTestCase):
  def setUp(self):
    # hg19-tRNA-Ala-TGC-3-1
    self.positions = annotate_positions('(.(.(.....(..........(.......(.(..............,....,..<..<<.<_............__........_____............._..._.._>.>.>>.........,.....<....<.<.<.<____.__............................................................................................_.>.>...>..>.>,..........,<<<<<<<_____>>>>>>...>,,..,<.<.....<.<.<._......._......._...._._...._._.............>.>>>....>............).)..)..).)........).)::::')

  def test_acceptor_with_gaps_on_both_sides(self):
    indices = [0, 2, 4, 10, 21, 29, 31]
    positions = ['1:397', '3:395', '5:386', '11:384', '22:381', '30:378', '32:376']
    sprinzl = ['1:72', '2:71', '3:70', '4:69', '5:68', '6:67', '7:66']
    self.testPositions(indices, positions, sprinzl)

  def test_dstem_with_gaps_on_both_sides(self):
    indices = [54, 57, 58, 60]
    positions = ['55:116', '58:115', '59:113', '61:111']
    sprinzl = ['10:25', '11:24', '12:23', '13:22']
    self.testPositions(indices, positions, sprinzl)

  def test_acstem_with_gaps_on_both_sides(self):
    indices = [127, 132, 134, 136, 138]
    positions = ['132:256', '137:254', '139:251', '141:247', '143:245']
    sprinzl = ['27:43', '28:42', '29:41', '30:40', '31:39']

class TestSSAlignmentAsp(AlignmentTestCase):
  def setUp(self):
    # hg19-tRNA-Asp-GTC-1-1
    self.positions = annotate_positions('(.(.(....(........(..(.(...........,..,..<<<.<_.________._..._>.>.>.>..,..<..<.<.<.<______................................_.>.>.....>..>.>,...,<<<<<<<_____>>>.>>>>,,..,<.<..<..<.<._.____._._.>.>>>....>...........).).)).)......).)::::')

  def test_acstem_with_gaps_on_both_sides(self):
    indices = [70, 73, 75, 77, 79]
    positions = ['75:138', '78:136', '80:133', '82:127', '84:125']
    sprinzl = ['27:43', '28:42', '29:41', '30:40', '31:39']
    self.testPositions(indices, positions, sprinzl)

if __name__ == '__main__':
    unittest.main()