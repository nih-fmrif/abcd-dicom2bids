def ftq_map(qc, directory):
    '''
    Given a list of ftq_series_ids from the abcd_fastqc01.txt file
    Return an expected BIDS directory and file outputs, as well as
    a mapping from ftq_series_id to resulting files
    
    Usage:
    outputs, mapping, errors, have_df, missing = ftq_map(qc_clean['ftq_series_id'], directory)
    
    Parameters
    ----------
    qc : list of strings
        list of ftq_series_id values
    directory : string
        rawdata directory to check expected files against
        
    Returns
    -------
    fin : dict
        dictionary of the resulting BIDS directory and file structure
        eg, fin[subject][session][anat] contains a list of expected files
    mapping : dict
        dictionary where the fqt_series_id is the key and the value is a list of files that should
        follow from the unpacking and converstion of the files in that series ID
    errors : list of strings
        list of of ftq_series_ids that did not parse according to len(str.split('_')) == 4
    have_df : data frame
        data frame consisting of the ftq_series_ids and whether the code was able to find the corresponding files
    missing : list of strings
        list of expected files the code could not find on disk
    '''
    fin = {} # will contain the expected directory and file structure
    mapping = {} # will contain the mapping of ftq_series_id to files
    errors = [] # list of ftq_series_ids that didn't parse 
    done_ftq = [] # list of ftq_ids
    done_yesno = [] # yes/no (1/0) for ftq_ids
    missing = [] # list of missing files
    
    # split ftq_series_id into constituent parts
    ftq_lst = []
    for i in qc:
        temp = i.split('_')
        if len(temp) == 4:
            temp.append(i)
            ftq_lst.append(temp)
        else:
            errors.append(i)
    ftq_df = pd.DataFrame(ftq_lst, columns=['sub', 'ses', 'series', 'time', 'ftq_series_id'])
    
    # subjects
    for s in ftq_df['sub'].unique():
        # s is subject ID w/o 'sub-' and sub has 'sub-'
        sub = f'sub-{s}'
        fin[sub] = {}
        # filter ftq_series_ids to grab this subject
        sub_df = ftq_df[ftq_df['sub'] == s]
        
        # sessions
        # ss is session w/o 'ses-' and ses has 'ses-'
        for ss in sub_df['ses'].unique():
            ses = f'ses-{ss}'
            fin[sub][ses] = {'anat':[], 'dwi': [], 'fmap': [], 'func': []}
            
            # filter ftq_series_ids to grab this session
            ses_df = sub_df[sub_df['ses'] == ss]
            
            # series
            for series in ses_df['series'].unique():
                # filter ftq_series_ids to grab this series
                run_df = ses_df[ses_df['series'] == series]
                # make sure then are sorted (for chronology)
                run_df = run_df.sort_values(by='time')
                
                # grab run index, all rows will have a run-xx assignment, but only > 1 runs will 
                # have 'run-' included in filenames
                runs = []
                for r in range(run_df.shape[0]):
                    runs.append(f'{r+1}'.zfill(2))
                run_df['run'] = runs
                
                if series == 'ABCD-T1':
                    # only 1 run
                    if run_df.shape[0] == 1:
                        # prep the ftq_series_id as a string, so it can be a key in the mapping
                        ftq = list(run_df.loc[:,'ftq_series_id'])[0]
                        mapping[ftq] = [(f'{sub}_{ses}_T1w.nii.gz')]
                        mapping[ftq].append(f'{sub}_{ses}_T1w.json')
                        fin[sub][ses]['anat'].extend(mapping[ftq])
                        # check if we have it
                        done_ftq.append(ftq)
                        have = 1
                        for f in mapping[ftq]:
                            if not os.path.isfile(f'{directory}/{sub}/{ses}/anat/{f}'):
                                have = 0
                                missing.append(f)
                        done_yesno.append(have)
                    # more than 1 run
                    else:
                        for i, r in enumerate(run_df['run']):
                            ftq = list(run_df.loc[:,'ftq_series_id'])[i]
                            mapping[ftq] = [(f'{sub}_{ses}_run-{r}_T1w.nii.gz')]
                            mapping[ftq].append(f'{sub}_{ses}_run-{r}_T1w.json')
                            fin[sub][ses]['anat'].extend(mapping[ftq])
                            # check if we have it
                            done_ftq.append(ftq)
                            have = 1
                            for f in mapping[ftq]:
                                if not os.path.isfile(f'{directory}/{sub}/{ses}/anat/{f}'):
                                    have = 0
                                    missing.append(f)
                            done_yesno.append(have)
                    
                elif series == 'ABCD-T1-NORM':
                    if run_df.shape[0] == 1:
                        ftq = list(run_df.loc[:,'ftq_series_id'])[0]
                        mapping[ftq] = [(f'{sub}_{ses}_rec-normalized_T1w.nii.gz')]
                        mapping[ftq].append(f'{sub}_{ses}_rec-normalized_T1w.json')
                        fin[sub][ses]['anat'].extend(mapping[ftq])
                        # check if we have it
                        done_ftq.append(ftq)
                        have = 1
                        for f in mapping[ftq]:
                            if not os.path.isfile(f'{directory}/{sub}/{ses}/anat/{f}'):
                                have = 0
                                missing.append(f)
                        done_yesno.append(have)
                    else:
                        for i, r in enumerate(run_df['run']):
                            ftq = list(run_df.loc[:,'ftq_series_id'])[i]
                            mapping[ftq] = [(f'{sub}_{ses}_rec-normalized_run-{r}_T1w.nii.gz')]
                            mapping[ftq].append(f'{sub}_{ses}_rec-normalized_run-{r}_T1w.json')
                            fin[sub][ses]['anat'].extend(mapping[ftq])
                            # check if we have it
                            done_ftq.append(ftq)
                            have = 1
                            for f in mapping[ftq]:
                                if not os.path.isfile(f'{directory}/{sub}/{ses}/anat/{f}'):
                                    have = 0
                                    missing.append(f)
                            done_yesno.append(have)
                    
                elif series == 'ABCD-T2':
                    if run_df.shape[0] == 1:
                        ftq = list(run_df.loc[:,'ftq_series_id'])[0]
                        mapping[ftq] = [(f'{sub}_{ses}_T2w.nii.gz')]
                        mapping[ftq].append(f'{sub}_{ses}_T2w.json')
                        fin[sub][ses]['anat'].extend(mapping[ftq])
                        # check if we have it
                        done_ftq.append(ftq)
                        have = 1
                        for f in mapping[ftq]:
                            if not os.path.isfile(f'{directory}/{sub}/{ses}/anat/{f}'):
                                have = 0
                                missing.append(f)
                        done_yesno.append(have)
                    else:
                        for i, r in enumerate(run_df['run']):
                            ftq = list(run_df.loc[:,'ftq_series_id'])[i]
                            mapping[ftq] = [(f'{sub}_{ses}_run-{r}_T2w.nii.gz')]
                            mapping[ftq].append(f'{sub}_{ses}_run-{r}_T2w.json')
                            fin[sub][ses]['anat'].extend(mapping[ftq])
                            # check if we have it
                            done_ftq.append(ftq)
                            have = 1
                            for f in mapping[ftq]:
                                if not os.path.isfile(f'{directory}/{sub}/{ses}/anat/{f}'):
                                    have = 0
                                    missing.append(f)
                            done_yesno.append(have)
                    
                elif series == 'ABCD-T2-NORM':
                    if run_df.shape[0] == 1:
                        ftq = list(run_df.loc[:,'ftq_series_id'])[0]
                        mapping[ftq] = [(f'{sub}_{ses}_rec-normalized_T2w.nii.gz')]
                        mapping[ftq].append(f'{sub}_{ses}_rec-normalized_T2w.json')
                        fin[sub][ses]['anat'].extend(mapping[ftq])
                        # check if we have it
                        done_ftq.append(ftq)
                        have = 1
                        for f in mapping[ftq]:
                            if not os.path.isfile(f'{directory}/{sub}/{ses}/anat/{f}'):
                                have = 0
                                missing.append(f)
                        done_yesno.append(have)
                    else:
                        for i, r in enumerate(run_df['run']):
                            ftq = list(run_df.loc[:,'ftq_series_id'])[i]
                            mapping[ftq] = [(f'{sub}_{ses}_rec-normalized_run-{r}_T2w.nii.gz')]
                            mapping[ftq].append(f'{sub}_{ses}_rec-normalized_run-{r}_T2w.json')
                            fin[sub][ses]['anat'].extend(mapping[ftq])
                            # check if we have it
                            done_ftq.append(ftq)
                            have = 1
                            for f in mapping[ftq]:
                                if not os.path.isfile(f'{directory}/{sub}/{ses}/anat/{f}'):
                                    have = 0
                                    missing.append(f)
                            done_yesno.append(have)
                    
                elif series == 'ABCD-DTI':
                    if run_df.shape[0] == 1:
                        ftq = list(run_df.loc[:,'ftq_series_id'])[0]
                        mapping[ftq] = [(f'{sub}_{ses}_dwi.nii.gz')]
                        mapping[ftq].append(f'{sub}_{ses}_dwi.json')
                        mapping[ftq].append(f'{sub}_{ses}_dwi.bval')
                        mapping[ftq].append(f'{sub}_{ses}_dwi.bvec')
                        fin[sub][ses]['dwi'].extend(mapping[ftq])
                        # check if we have it
                        done_ftq.append(ftq)
                        have = 1
                        for f in mapping[ftq]:
                            if not os.path.isfile(f'{directory}/{sub}/{ses}/dwi/{f}'):
                                have = 0
                                missing.append(f)
                        done_yesno.append(have)
                    else:
                        for i, r in enumerate(run_df['run']):
                            ftq = list(run_df.loc[:,'ftq_series_id'])[i]
                            mapping[ftq] = [(f'{sub}_{ses}_run-{r}_dwi.nii.gz')]
                            mapping[ftq].append(f'{sub}_{ses}_run-{r}_dwi.json')
                            mapping[ftq].append(f'{sub}_{ses}_run-{r}_dwi.bval')
                            mapping[ftq].append(f'{sub}_{ses}_run-{r}_dwi.bvec')
                            fin[sub][ses]['dwi'].extend(mapping[ftq])
                            # check if we have it
                            done_ftq.append(ftq)
                            have = 1
                            for f in mapping[ftq]:
                                if not os.path.isfile(f'{directory}/{sub}/{ses}/dwi/{f}'):
                                    have = 0
                                    missing.append(f)
                            done_yesno.append(have)
                    
                elif series == 'ABCD-MID-fMRI':
                    if run_df.shape[0] == 1:
                        ftq = list(run_df.loc[:,'ftq_series_id'])[0]
                        mapping[ftq] = [(f'{sub}_{ses}_task-MID_bold.nii.gz')]
                        mapping[ftq].append(f'{sub}_{ses}_task-MID_bold.json')
                        fin[sub][ses]['func'].extend(mapping[ftq])
                        # check if we have it
                        done_ftq.append(ftq)
                        have = 1
                        for f in mapping[ftq]:
                            if not os.path.isfile(f'{directory}/{sub}/{ses}/func/{f}'):
                                have = 0
                                missing.append(f)
                        done_yesno.append(have)
                    else:
                        for i, r in enumerate(run_df['run']):
                            ftq = list(run_df.loc[:,'ftq_series_id'])[i]
                            mapping[ftq] = [(f'{sub}_{ses}_task-MID_run-{r}_bold.nii.gz')]
                            mapping[ftq].append(f'{sub}_{ses}_task-MID_run-{r}_bold.json')
                            fin[sub][ses]['func'].extend(mapping[ftq])
                            # check if we have it
                            done_ftq.append(ftq)
                            have = 1
                            for f in mapping[ftq]:
                                if not os.path.isfile(f'{directory}/{sub}/{ses}/func/{f}'):
                                    have = 0
                                    missing.append(f)
                            done_yesno.append(have)
                            
                elif series == 'ABCD-nBack-fMRI':
                    if run_df.shape[0] == 1:
                        ftq = list(run_df.loc[:,'ftq_series_id'])[0]
                        mapping[ftq] = [(f'{sub}_{ses}_task-nback_bold.nii.gz')]
                        mapping[ftq].append(f'{sub}_{ses}_task-nback_bold.json')
                        fin[sub][ses]['func'].extend(mapping[ftq])
                        # check if we have it
                        done_ftq.append(ftq)
                        have = 1
                        for f in mapping[ftq]:
                            if not os.path.isfile(f'{directory}/{sub}/{ses}/func/{f}'):
                                have = 0
                                missing.append(f)
                        done_yesno.append(have)
                    else:
                        for i, r in enumerate(run_df['run']):
                            ftq = list(run_df.loc[:,'ftq_series_id'])[i]
                            mapping[ftq] = [(f'{sub}_{ses}_task-nback_run-{r}_bold.nii.gz')]
                            mapping[ftq].append(f'{sub}_{ses}_task-nback_run-{r}_bold.json')
                            fin[sub][ses]['func'].extend(mapping[ftq])
                            # check if we have it
                            done_ftq.append(ftq)
                            have = 1
                            for f in mapping[ftq]:
                                if not os.path.isfile(f'{directory}/{sub}/{ses}/func/{f}'):
                                    have = 0
                                    missing.append(f)
                            done_yesno.append(have)
                            
                elif series == 'ABCD-SST-fMRI':
                    if run_df.shape[0] == 1:
                        ftq = list(run_df.loc[:,'ftq_series_id'])[0]
                        mapping[ftq] = [(f'{sub}_{ses}_task-SST_bold.nii.gz')]
                        mapping[ftq].append(f'{sub}_{ses}_task-SST_bold.json')
                        fin[sub][ses]['func'].extend(mapping[ftq])
                        # check if we have it
                        done_ftq.append(ftq)
                        have = 1
                        for f in mapping[ftq]:
                            if not os.path.isfile(f'{directory}/{sub}/{ses}/func/{f}'):
                                have = 0
                                missing.append(f)
                        done_yesno.append(have)
                    else:
                        for i, r in enumerate(run_df['run']):
                            ftq = list(run_df.loc[:,'ftq_series_id'])[i]
                            mapping[ftq] = [(f'{sub}_{ses}_task-SST_run-{r}_bold.nii.gz')]
                            mapping[ftq].append(f'{sub}_{ses}_task-SST_run-{r}_bold.json')
                            fin[sub][ses]['func'].extend(mapping[ftq])
                            # check if we have it
                            done_ftq.append(ftq)
                            have = 1
                            for f in mapping[ftq]:
                                if not os.path.isfile(f'{directory}/{sub}/{ses}/func/{f}'):
                                    have = 0
                                    missing.append(f)
                            done_yesno.append(have)
                            
                elif series == 'ABCD-rsfMRI':
                    if run_df.shape[0] == 1:
                        ftq = list(run_df.loc[:,'ftq_series_id'])[0]
                        mapping[ftq] = [(f'{sub}_{ses}_task-rest_bold.nii.gz')]
                        mapping[ftq].append(f'{sub}_{ses}_task-rest_bold.json')
                        fin[sub][ses]['func'].extend(mapping[ftq])
                        # check if we have it
                        done_ftq.append(ftq)
                        have = 1
                        for f in mapping[ftq]:
                            if not os.path.isfile(f'{directory}/{sub}/{ses}/func/{f}'):
                                have = 0
                                missing.append(f)
                        done_yesno.append(have)
                    else:
                        for i, r in enumerate(run_df['run']):
                            ftq = list(run_df.loc[:,'ftq_series_id'])[i]
                            mapping[ftq] = [(f'{sub}_{ses}_task-rest_run-{r}_bold.nii.gz')]
                            mapping[ftq].append(f'{sub}_{ses}_task-rest_run-{r}_bold.json')
                            fin[sub][ses]['func'].extend(mapping[ftq])
                            # check if we have it
                            done_ftq.append(ftq)
                            have = 1
                            for f in mapping[ftq]:
                                if not os.path.isfile(f'{directory}/{sub}/{ses}/func/{f}'):
                                    have = 0
                                    missing.append(f)
                            done_yesno.append(have)
                                                      
    have_df = pd.DataFrame(list(zip(done_ftq, done_yesno)), columns=['ftq', 'done'])
                            
    return fin, mapping, errors, have_df, missing