import md
import os.path

md.run_md()

if os.path.isfile('cu.traj'):
    print ("File exist")
else:
    print ("File not exist")
