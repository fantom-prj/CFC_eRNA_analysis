import sys
import os
os.environ["CUDA_VISIBLE_DEVICES"] = sys.argv[1]
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
import numpy as np
import pandas as pd

# python3 junction_run.py 0 all_ontCAGE.SS3J41.tsv >all_ontCAGE.SS3J41.tsv.out 2>&1 &
context = 10000
paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
models = [load_model(resource_filename('spliceai', x)) for x in paths]

columns = ['ID', 'acceptor max', 'donor max', 'acceptor', 'donor']
junction_folder = "junction"
output_folder = "spliceai_out"
# for filename in os.listdir(junction_folder):
#     if filename.endswith(".tsv"):
filename = sys.argv[2]
print("Starting " + filename)
new_rows = []
info = filename.split(".")
info = info[-2]
info = info.split("J")
ss = info[0]
d = int(info[1])
df = pd.read_csv(os.path.join(junction_folder, filename), sep='\t', header=None, names=['id', 'seq'])
for index, row in df.iterrows():
    if index % 100 == 0:
        print(index, end=" ", flush=True)
    x = one_hot_encode('N'*(context//2) + row['seq'] + 'N'*(context//2))[None, :]
    y = np.mean([models[m].predict(x) for m in range(5)], axis=0)
    # y = models[-1].predict(x)
    acceptor_prob = y[0, :, 1]
    donor_prob = y[0, :, 2]
    new_row = [row['id'], np.amax(acceptor_prob), np.amax(donor_prob), ','.join(acceptor_prob.astype(str)), ','.join(donor_prob.astype(str))]
    new_rows.append(new_row)

df = pd.DataFrame(new_rows, columns=columns)
df.to_csv(output_folder + '/' + os.path.splitext(filename)[0] + "_5SpliceAI.tsv", sep='\t', index=False)
print("\nFinished "+ filename)

