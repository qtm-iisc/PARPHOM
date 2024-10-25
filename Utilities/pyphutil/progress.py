"""
    Written by Shinjan Mandal as an utility for PARPHOM
    Github repo: https://github.com/qtm-iisc/PARPHOM
    Email : mshinjan@iisc.ac.in
"""
import sys

def printProgressBar(i,m,postText,n_bar=20):
    """
        Prints progress of the calculations
    """
    j=i/m
    sys.stdout.write('\r')
    sys.stdout.write(f" {'█' * int(n_bar * j):{n_bar}s} {int(100 * j)}% {postText}")
    sys.stdout.flush()
    return
