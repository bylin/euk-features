import sys
import re

scores = {}
positions = {'2': '73', '3': '1:72', '4': '2:71', '5': '3:70', '6': '4:69', '7': '5:68', '8': '6:67', '9': '7:66', '10': '8', '11': '9', '16': '10:25', '17': '11:24', '18': '12:23', '19': '13:22', '20': '14', '21': '15', '22': '16', '23': '17', '24': '17a', '25': '18', '26': '19', '27': '20', '28': '20a', '29': '20b', '30': '21', '33': '26', '34': '27:43', '35': '28:42', '36': '29:41', '37': '30:40', '38': '31:39', '39': '32', '40': '33', '41': '34', '42': '35', '43': '36', '44': '37', '45': '38', '48': '44', '49': '45', '52': 'V11:V21', '53': 'V12:V22', '54': 'V13:V23', '55': 'V14:V24', '56': 'V15:V25', '57': 'V16:V26', '58': 'V17:V27', '59': 'V1', '60': 'V2', '61': 'V3', '62': 'V4', '63': 'V5', '66': '46', '67': '47', '68': '48', '69': '49:65', '70': '50:64', '71': '51:63', '72': '52:62', '73': '53:61', '74': '54', '75': '55', '76': '56', '77': '57', '78': '58', '79': '59', '80': '60'}
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
  if rowid == "81":
    break

  # skip special node rows
  if int(rowid) in [0, 1, 12, 13, 14, 15, 31, 32, 46, 47, 50, 51, 64, 65]:
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

for score in sorted(scores): print('{}\t{}'.format(score, scores[score]))