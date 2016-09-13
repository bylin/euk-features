import re
from collections import Counter
from math import log2

def get_positions(input_file):
  alignment_fhandle = open(input_file)
  positions = [] # list containing each position in the tRNA
  # first, get secondary structure
  for line in alignment_fhandle:
    if line[0:12] == '#=GC SS_cons':
      ss = line.strip().split()[-1]
  # parse secondary structure into regions and positions
  positions = annotate_positions(ss)
  # get counts for each position by parsing Stockholm file
  positions = count_positions(input_file, positions)
  return positions

class Position:
  def __init__(self, position, region, region_number, index=-1, paired=False, counts=None, num_obs=0, entropy=-1):
    self.position = position
    self.region = region
    self.region_number = region_number
    self.index = index
    self.paired = paired
    if counts is None:
      self.counts = Counter()
    else:
      self.counts = counts
    self.num_obs = num_obs
    self.entropy = entropy
  
  def __str__(self):
    return "Position {} ({} #{})".format(self.position, self.region, self.region_number)

class Region:
  def __init__(self, lower, upper, name):
    self.lower = lower
    self.upper = upper
    self.name = name
  
  def __str__(self):
    return "Region {} ({} ~ {})".format(self.name, self.lower, self.upper)

def annotate_positions(ss):
  loop_indices = [r.span() for r in re.finditer('\(+|<+|_+|>+|\)+', ss)]
  if len(loop_indices) == 14: regions = ['acceptor', 'dstem', 'dloop', 'dstem', 'acstem', 'acloop', 'acstem', 'vstem', 'vloop', 'vstem', 'tpcstem', 'tpcloop', 'tpcstem', 'acstem']
  elif len(loop_indices) == 11: regions = ['acceptor', 'dstem', 'dloop', 'dstem', 'acstem', 'acloop', 'acstem', 'tpcstem', 'tpcloop', 'tpcstem', 'acstem']
  regions = [Region(indices[0], indices[1], name) for indices, name in zip(loop_indices, regions[:])]
  region = regions[0]
  positions = []
  region_index = 0 # index to be used to iterate through regions list
  region_numbering = 0 # base numbering within a region
  for position in range(len(ss)):
    if region_index < len(regions):
      region = regions[region_index]
    else:
      positions.append(Position(position=str(position + 1), region='single', region_number=region_numbering + 1, index=len(positions), paired=False))
      continue
    if position < region.lower: # before the next region starts (or if it's the last region), annotate as single bases
      positions.append(Position(position=str(position + 1), region='single', region_number=region_numbering + 1, index=len(positions), paired=False))
      region_numbering += 1
    elif position == region.lower: # start of region: begin region numbering at 1
      if ss[position] == "(":
        paired_base = regions[-1].upper - 1
        positions.append(Position(position='{}:{}'.format(position + 1, paired_base + 1), region=region.name, region_number=1, index=len(positions), paired=True))
      elif ss[position] == "<":
        paired_base = regions[region_index + 2].upper - 1
        positions.append(Position(position='{}:{}'.format(position + 1, paired_base + 1), region=region.name, region_number=1, index=len(positions), paired=True))
      elif ss[position] in [')', '>']:
        pass
      else:
        positions.append(Position(position=str(position + 1), region=region.name, region_number=region_numbering + 1, index=len(positions), paired=False))
      region_numbering = 1
    elif position > region.lower and position <= region.upper - 1: # inside region: increment region numbering normally
      # find paired base, or skip if base is the opposite strand
      if ss[position] == "(":
        paired_base = regions[-1].upper - (position - regions[0].lower) - 1
        positions.append(Position(position='{}:{}'.format(position + 1, paired_base + 1), region=region.name, region_number=region_numbering + 1, index=len(positions), paired=True))
      elif ss[position] == "<":
        paired_base = regions[region_index + 2].upper - region_numbering - 1
        positions.append(Position(position='{}:{}'.format(position + 1, paired_base + 1), region=region.name, region_number=region_numbering + 1, index=len(positions), paired=True))
      elif ss[position] in [')', '>']:
        pass
      else:
        positions.append(Position(position=str(position + 1), region=region.name, region_number=region_numbering+1, index=len(positions), paired=False))
      region_numbering += 1
    if position == region.upper - 1: # end of region, reset region index and increment region number
      region_index += 1
      region_numbering = 0
    
  return positions

def count_positions(input_file, positions):
  # We've now annotated all of the positions. Time to count bases at each position.
  # Loop through all sequences in alignment
  alignment_fhandle = open(input_file) # refresh handle
  for line in alignment_fhandle:
    if line[0] in ["#", '\n', '/']: continue
    if len(line.split()) == 1: print(line)
    seqname, seq = line.strip().split()
    for position_index, position in enumerate(positions):
      if position.paired:
        index1, index2 = position.position.split(':')
        index1, index2 = int(index1), int(index2)
        base_pair = "{}:{}".format(seq[index1 - 1], seq[index2 - 1])
        positions[position_index].counts[base_pair] += 1
      else:
        index = int(position.position)
        base = seq[index - 1]
        positions[position_index].counts[base] += 1
  for index, position in enumerate(positions):
    positions[index].num_obs = sum(position.counts.values())
    positions[index].entropy = calculate_positional_entropy(positions[index])
  return positions

def calculate_positional_entropy(position):
  return sum(-count / position.num_obs * log2(count / position.num_obs) for count in position.counts.values())

def position_generator(positions, threshold=0.98):
  for position in positions: # iterate through positions
    if position.paired:
      combos = [{'A:U': ['A:U'], 'U:A': ['U:A'], 'G:C': ['G:C'], 'C:G': ['C:G'], 'G:U': ['G:U'], 'U:G': ['U:G'], 'A:A': ['A:A'], 'C:C': ['C:C'], 'G:G': ['G:G'], 'U:U': ['U:U'], 'A:G': ['A:G'], 'G:A': ['G:A'], 'A:C': ['A:C'], 'C:A': ['C:A'], 'C:U': ['C:U'], 'U:C': ['U:C']},
                {'R:Y': ['G:C', 'A:U'], 'Y:R': ['C:G', 'U:A'], 'S:S': ['G:C', 'C:G'], 'W:W': ['A:U', 'U:A'], 'W:O': ['G:U', 'U:G']},
                {'B:V': ['G:C', 'C:G', 'U:A'], 'V:B': ['G:C', 'C:G', 'A:U'], 'D:H': ['A:U', 'U:A', 'G:C'], 'H:D': ['A:U', 'U:A', 'C:G']},
                {'W:C': ['A:U', 'U:A', 'G:C', 'C:G']},
                {'G:N': ['G:C', 'G:U', 'C:G', 'U:G'], 'U:N': ['U:G', 'U:A', 'A:U', 'G:U']},
                {'C:O': ['A:U', 'U:A', 'G:C', 'C:G', 'G:U', 'U:G']},
                {'Y:Y': ['U:C', 'C:U', 'U:U', 'C:C'], 'R:R': ['A:G', 'G:A', 'A:A', 'G:G']},
                {'M:M': ['A:A', 'G:G', 'C:C', 'U:U', 'A:G', 'A:C', 'C:A', 'C:U', 'G:A', 'U:C']}]
    else: 
      combos = [{'A': ['A'], 'C': ['C'], 'G': ['G'], 'U': ['U']},
                {'R': ['A', 'G'], 'Y': ['C', 'U'], 'K': ['G', 'U'], 'M': ['A', 'C'], 'S': ['G', 'C'], 'W': ['A', 'U']},
                {'B': ['C', 'G', 'U'], 'D': ['G', 'A', 'U'], 'H': ['A', 'C', 'U'], 'V': ['G', 'C', 'A']}]
      # Single base combinations. Three dictionaries, organized by rank.
    max_freq = 0
    best_symbol = ''
    for rank_dict in combos:
      # Each dictionary contains symbols (keys) and possible bases or base pairs the symbol resolves to (values).
      # Check which symbol has the highest frequency.
      for symbol, bases in rank_dict.items():
        freqs = [position.counts[base] for base in bases] / sum(position.counts.values())
        # Each possible base/bp needs to occur in at least 2% of tRNAs
        if any(freqs) < 0.02: continue
        freq = sum(freqs)
        if max_freq < freq:
          max_freq = freq
          best_symbol = symbol
      if max_freq > threshold:
        yield position, best_symbol, max_freq
        break
    if max_freq < threshold:
      yield position, best_symbol, max_freq
