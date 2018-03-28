"""
Download all of the observational data we need
"""
import os
import drms

EMAIL = 'wtb2@rice.edu'
DOWNLOAD_DIR = '/storage-home/w/wtb2/data/timelag_synthesis_v2/observational_data/aia'

def download_data(channel, index=None):
    """
    Download AR 1158 AIA data from JSOC using drms client
    """
    c = drms.Client(email=EMAIL, verbose=True)
    # Time over which to iterate
    start_time = '2011-02-12T9:32:13' # -interval/2 from when the AR was observed
    interval = '12h'
    cadence = '12s'
    # check that download dir is available
    if not os.path.exists(DOWNLOAD_DIR):
        os.mkdir(DOWNLOAD_DIR)
    # build query
    q = f'aia.lev1_euv_12s[{start_time}/{interval}@{cadence}][{channel}]{{image}}'
    r = c.export(q, method='url_quick')
    if r.status != 0:
        raise ValueError(f'Bad query with exit status {r.status}')
    # download
    df = r.download(DOWNLOAD_DIR, index=index, fname_from_rec=True)

    return df

if __name__ == '__main__':
    channels = [94,131,171,193,211,335]
    for chan in channels:
        df = download_data(chan)
        print(df)
