#! python3
# batch_bader_analysis.py - do bader analysis from all subfolders in the current working folder.
# usage: python3 batch_bader_analysis.py
# author - Jiangcg
import os
from alive_progress import alive_bar

count = 0
global current_path
current_path = os.getcwd()
for foldername, subfolders, filenames in os.walk(current_path):
    with alive_bar(len(subfolders), title='Writing to DATABASE') as bar:
        for i in range(len(subfolders)):
            os.chdir(subfolders[i])
            os.system('dobader')
            bar()
            count += 1
            print('Bader electron analysis of ' + subfolders[i] + ' has been done! ')
            os.chdir(current_path)
        break
print( str(count) + ' folders bader elctron analysis have been done! Good luck!' )
