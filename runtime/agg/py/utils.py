import sys, subprocess, os

EOF_IS_NL = True

def read_file(fname):
  try:
    return open(fname, 'r').readlines()
  except IOError as _err:
    # os.path.basename(sys.argv[0]) + ": " + 
    sys.stderr.write(f + ": " + _err.strerror + "\n")

def read_all(): 
  global EOF_IS_NL
  all_contents = []
  for f in sys.argv[1:]:
    contents = read_file(f)
    EOF_IS_NL = EOF_IS_NL and contents[-1].endswith('\n')
    all_contents.append(contents)
  return all_contents

#def out():
#  if (sys.version_info.major == 2):
#    print **locals(),
#  else:
#    print(**locals(), end=' ')
def out(s):
  global EOF_IS_NL
  if not s.endswith('\n') and EOF_IS_NL:
    sys.stdout.write(s + '\n')
  else:
    sys.stdout.write(s)
  sys.stdout.flush()

def cmd():
  c = sys.argv[0].replace(".py", "").replace("./", "")
  if 'MAP_CMD' in os.environ:
    c = os.environ['MAP_CMD']
  return c

def execute(command, data):
    p = subprocess.Popen([command], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    return p.communicate(data)[0]
    # Python 3 equivalent:
    # p = subprocess.run([cmd], stdout=subprocess.PIPE, input=data, encoding='ascii', check=True)
    # return p.stdout

def help(c=cmd()):
  m = c
  s = ' <(cat 1.t | ' + m + ') <(cat 2.t | ' + m + ') <(cat 3.t | ' + m + ')'
  if len(sys.argv) < 2:
    print 'echo "one\\ntwo\\nthree" > 1.t'
    print 'echo "four\\nfive" > 2.t'
    print 'echo "one\\ntwo\\nthree" > 3.t'
    print sys.argv[0] + s
    print '\nTo quickly test equivalence: '
    print sys.argv[0] + s + ' > par.txt'
    print 'cat 1.t 2.t 3.t | ' + m + ' > seq.txt'
    print 'diff {seq,par}.txt'


