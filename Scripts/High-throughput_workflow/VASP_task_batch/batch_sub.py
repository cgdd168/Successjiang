import os

current_path = os.getcwd()

for root, dirs, files in os.walk(current_path):
    if 'script' in files:
        num = os.path.join(root, 'script')
        os.chdir(root)
        os.system('qsub script')