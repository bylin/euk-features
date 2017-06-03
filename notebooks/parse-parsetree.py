''' Helper script for analyzing tRNA sequence features. Parses single tRNA .tfile parsetrees. '''
import sys
import re

scores = {}
identities = {}
positions = {4: '73', 5: '1:72', 6: '2:71', 7: '3:70', 8: '4:69', 9: '5:68', 10: '6:67', 11: '7:66', 12: '8', 13: '9', 18: '10:25', 19: '11:24', 20: '12:23', 21: '13:22', 22: '14', 23: '15', 24: '16', 25: '17', 26: '17a', 27: '18', 28: '19', 29: '20', 30: '20a', 31: '20b', 32: '21', 35: '26', 36: '27:43', 37: '28:42', 38: '29:41', 39: '30:40', 40: '31:39', 41: '32', 42: '33', 43: '34', 44: '35', 45: '36', 46: '37', 47: '38', 50: '44', 51: '45', 54: 'V11:V21', 55: 'V12:V22', 56: 'V13:V23', 57: 'V14:V24', 58: 'V15:V25', 59: 'V16:V26', 60: 'V17:V27', 61: 'V1', 62: 'V2', 63: 'V3', 64: 'V4', 65: 'V5', 68: '46', 69: '47', 70: '48', 71: '49:65', 72: '50:64', 73: '51:63', 74: '52:62', 75: '53:61', 76: '54', 77: '55', 78: '56', 79: '57', 80: '58', 81: '59', 82: '60'}
skip_positions = [0, 1, 2, 3, 14, 15, 16, 17, 33, 34, 48, 49, 52, 53, 66, 67]
terminal_position = 83
doneParsingHeader = False

for line in open(sys.argv[1]):
  cols = line.strip().split()
  if len(cols) > 0 and cols[0] == '0':
    doneParsingHeader = True
    tsc = float(cols[8])
    continue
  if not doneParsingHeader:
    continue

  # parse row. columns: rowid, emitl, emitr, state, mode, nxtl, nxtr, prv, tsc, esc
  rowid = int(cols[0])
  emitl = cols[1]
  emitr = cols[2]
  state = cols[3]
  prev_tsc = tsc
  tsc = float(cols[8])
  esc = float(cols[9])

  # exit on terminal node
  if rowid == terminal_position:
    break

  # skip special node rows
  if rowid in skip_positions:
    tsc = float(cols[8])
    continue

  # add standard match positions to scores dict
  if state[-2:] in ["MR", "ML", "MP"]:
    scores[positions[rowid]] = prev_tsc + esc
    if state[-2:] == "MR": identities[positions[rowid]] = re.findall('[A-Z]', emitr)[0]
    if state[-2:] == "ML": identities[positions[rowid]] = re.findall('[A-Z]', emitl)[0]
    if state[-2:] == "MP": identities[positions[rowid]] = '{}:{}'.format(re.findall('[A-Z]', emitl)[0], re.findall('[A-Z]', emitr)[0])

  # for deletions, don't add the esc value
  if state[-1] == "D":
    scores[positions[rowid]] = prev_tsc
    identities[positions[rowid]] = '-'
    if ':' in positions[rowid]: identities[positions[rowid]] = '-:-'
    
  # for insertions, increment remaining positions by 1 and skip
  if state[-2:] in ["IL", "IR"]:
    terminal_position += 1
    for position in sorted(positions, reverse = True):
      if position < rowid: break
      positions[position + 1] = positions.pop(position)

    for i, position in reversed(list(enumerate(sorted(skip_positions)))):
      if position < rowid: break
      skip_positions[i] += 1



for position in sorted(scores): print('{}\t{}\t{}'.format(position, round(scores[position], 2), identities[position]))
