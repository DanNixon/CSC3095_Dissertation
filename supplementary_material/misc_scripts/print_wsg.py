# Prints the run title of all workspaces in a workspace group

wsg = 'MultiFiles'

for w in mtd[wsg]:
    print '{0} = {1}'.format(w, w.getRun().getProperty('run_title').valueAsStr)
