import sys
import re

scores = {}
positions = {'4': '73', '5': '1:72', '6': '2:71', '7': '3:70', '8': '4:69', '9': '5:68', '10': '6:67', '11': '7:66', '12': '8', '13': '9', '18': '10:25', '19': '11:24', '20': '12:23', '21': '13:22', '22': '14', '23': '15', '24': '16', '25': '17', '26': '17a', '27': '18', '28': '19', '29': '20', '30': '20a', '31': '20b', '32': '21', '35': '26', '36': '27:43', '37': '28:42', '38': '29:41', '39': '30:40', '40': '31:39', '41': '32', '42': '33', '43': '34', '44': '35', '45': '36', '46': '37', '47': '38', '50': '44', '51': '45', '54': 'V11:V21', '55': 'V12:V22', '56': 'V13:V23', '57': 'V14:V24', '58': 'V15:V25', '59': 'V16:V26', '60': 'V17:V27', '61': 'V1', '62': 'V2', '63': 'V3', '64': 'V4', '65': 'V5', '68': '46', '69': '47', '70': '48', '71': '49:65', '72': '50:64', '73': '51:63', '74': '52:62', '75': '53:61', '76': '54', '77': '55', '78': '56', '79': '57', '80': '58', '81': '59', '82': '60'}
doneParsingHeader = False

for line in open(sys.argv[1]):
  cols = line.strip().split()
  if len(cols) > 0 and cols[0] == '0':
    doneParsingHeader = True
  if not doneParsingHeader:
    continue

  # columns: rowid, emitl, emitr, state, mode, nxtl, nxtr, prv, tsc, esc

  rowid = cols[0]
  # exit on terminal node
  if rowid == "83":
    break

  # skip special node rows
  if int(rowid) in [0, 1, 2, 3, 14, 15, 16, 17, 33, 34, 48, 49, 52, 53, 66, 67]:
    tsc = float(cols[8])
    continue

  # parse each row
  state = cols[3]
  prev_tsc = tsc
  tsc = float(cols[8])
  esc = float(cols[9])
  if state[-2:] in ["MR", "ML", "MP"]:
    scores[positions[rowid]] = prev_tsc + esc

  # for deletions, don't add the tsc value
  if state[-1] == "D":
    tsc = float(cols[8])
    scores[positions[rowid]] = esc

for score in sorted(scores): print('{}\t{}'.format(score, round(scores[score], 2)))