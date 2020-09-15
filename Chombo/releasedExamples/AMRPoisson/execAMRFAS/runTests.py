# Create inputs files a series of AMR runs to test the error convergence
import re
import os

# Load up an inputs file and parse it into a dictionary
def readInputs(inputs_file):

    params = {}

    # Read in the file
    filedata = None
    with open(inputs_file, 'r') as f:
        filedata = f.readlines()

        for line in filedata:
            # print(line)
            # Ignore lines that are commented out
            if line.startswith('#'):
                continue

            # Remove anything after a #
            line = re.sub('#[^\n]*[\n]', '', line)
            # print(line)

            parts = line.split('=')
            if len(parts) > 1:
                key = parts[0].strip()
                val = parts[1].strip()
                params[key] = val

    # print(params)
    return params

def writeParamsFile(location, params, ignoreList = [], doSort=True):
    output_file = ''

    keyList = params.keys()
    if doSort:
        keyList.sort()    

        for key in keyList:
            if key not in ignoreList:
                output_file = output_file + '\n' + \
                key + '=' + str(params[key])


                with open(location, 'w') as f:
                    f.write(output_file)

def getRunTime(output_folder):
    timefile = os.path.join(output_folder, 'time.table.0')

   

    if os.path.exists(timefile):
       
        with open(timefile, 'r') as f:

            for line in f.readlines():
                
                match = re.findall('\[0\]main (\d+\.\d+)', line)
                if len(match) > 0:
                    runtime = float(match[0])
                    return runtime

            
            return float('NaN')

    else:
        print('Cannot open %s' % timefile)



# Start script
execFile = 'FASsolve2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex'
nxs = [16, 32, 64, 128, 256, 512]
params = readInputs('./inputs')
err = {}

output_text = '==================================\n'\
'Convergence test\n'\
'----------------------------------\n'\
'Nx   | Max err  | L1       | L2       | Rate  || Runtime (s)  | Rate \n'

for nx in nxs:

    params['grids.num_cells'] = str(nx) + " " + str(nx) + " " + str(nx)

    # Create grids
    params['grids.max_level'] = '2'
    params['grids.read_in_grids']       = 'true'
    params['grids.level_1_box_count']   = '1'
    bl = nx/4
    tr = bl + nx - 1
    params['grids.level_1_box_0_lo']   = str(bl) + " " + str(bl) + " " + str(bl)
    params['grids.level_1_box_0_hi']   = str(tr) + " " + str(tr) + " " + str(tr)

    bl = nx/2+nx/4
    tr = bl + nx-1
    params['grids.level_2_box_count']   = '1'
    params['grids.level_2_box_0_lo']  = str(bl) + " " + str(bl) + " " + str(bl)
    params['grids.level_2_box_0_hi']   = str(tr) + " " + str(tr) + " " + str(tr)


    # Print params file
    output_dir = './convergenceTest/' + str(nx) + '/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

        writeParamsFile(output_dir + 'inputs', params)

    # Do run
    cmd = "cd " + output_dir + "; " + "../../" + execFile + " inputs"
    os.system(cmd)

    # Collate results
    poutFile = output_dir + 'pout.0'

    runTime = getRunTime(output_dir)

    with open(poutFile, 'r') as f:
        filedata = f.readlines()

        for line in filedata:
            if not line[0:5] == "Error":
                continue

            # Get error
            floatRegex = '\d+\.\d+e-?\d+'
            regex = 'Error = ('+ floatRegex + ') \(max\), ('+ floatRegex + ') \(L1\), ('+ floatRegex + ') \(L2\)'

            results = re.findall(regex, line)

            if len(results) > 0:
                res = [float(r) for r in results[0]] 
                res.append(runTime)
                err[str(nx)] = res



                # Compute rate of convergence
                if nx > nxs[0]:

                    rate = (err[str(nx/2)][0]/err[str(nx)][0])/2
                    runTimeRate = (err[str(nx)][3]/err[str(nx/2)][3])/4
                else:
                    rate = 0.0
                    runTimeRate = 0.0

                output_line = '%-4d | %.2e | %.2e | %.2e | %.3f || %7.3f      | %1.3f \n' % (nx, res[0], res[1], res[2], rate, runTime, runTimeRate )
                output_text += output_line

print(output_text)
print('==================================')