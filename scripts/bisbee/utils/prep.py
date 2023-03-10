import pandas as pd
import sys, os, re
import h5py
import numpy as np

events_file=sys.argv[1]
event_type=sys.argv[2]
outname=sys.argv[3]
spladder_ver=int(sys.argv[4])
if len(sys.argv)>5:
    sample_file=sys.argv[5]
if len(sys.argv)>6:
    chunk_num=int(sys.argv[6])
    chunk_size=int(sys.argv[7])

event_dict={'alt_3prime':'alt_3prime','alt_5prime':'alt_5prime','exon_skip':'exon_skip','mult_exon_skip':'mult_exon_skip','intron_retention':'intron_retention','mutex_exons':'mutex_exons'}

f=h5py.File(events_file,'r')
try:
    start_pos=chunk_num*chunk_size
    end_pos=min((chunk_num+1)*chunk_size,len(f['gene_idx']))
except:
    start_pos=0
    end_pos=len(f['gene_idx'])

gene_idx=np.array(f['gene_idx'][start_pos:end_pos]).astype(int)
conf=f['confirmed'][start_pos:end_pos]
gene_names=np.array(list(map(lambda x: x.decode('UTF-8'),f['gene_names'])))
gene_chr=np.array(list(map(lambda x: x.decode('UTF-8'),f['gene_chr'])))
gene_strand=np.array(list(map(lambda x: x.decode('UTF-8'),f['gene_strand'])))
sample_names=np.array(list(map(lambda x: x.decode('UTF-8'),f['strains'])))
try:
    select_samples=pd.read_table(sample_file,usecols=[0],names=["samples"])
except:
    print('no sample file')
    select_samples=pd.DataFrame(sample_names,columns=["samples"])
samples,idx1,idx2=np.intersect1d(sample_names,select_samples.samples,return_indices=True)
idx1.sort()


data=pd.DataFrame({'gene':gene_names[gene_idx],'strand':gene_strand[gene_idx],'contig':gene_chr[gene_idx]})
data['event_id']=list(map(lambda x: event_dict[event_type] + '_' + str(x),range(start_pos+1,end_pos+1)))
data['confirmed']=conf

if event_type=='mutex_exons':
    event_pos=f['event_pos'][start_pos:end_pos,:,:]
else:
    event_pos=f['event_pos'][start_pos:end_pos,:]

if event_type=='alt_3prime':
    data['exon_const_start']=event_pos[:,0]+1
    data['exon_const_end']=event_pos[:,1]
    data['exon_alt1_start']=event_pos[:,4]+1
    data['exon_alt1_end']=event_pos[:,5]
    data['exon_alt2_start']=event_pos[:,6]+1
    data['exon_alt2_end']=event_pos[:,7]
    iso1_str=data.apply(lambda x: str(x['exon_const_end'])+'j'+str(x['exon_alt2_start']),axis=1)
    iso2_str=data.apply(lambda x: str(x['exon_const_end'])+'j'+str(x['exon_alt1_start']),axis=1)
    if sum(data['strand']=='-')>0:
        data.loc[data['strand']=='-','exon_const_start']=event_pos[data['strand']=='-',4]+1
        data.loc[data['strand']=='-','exon_const_end']=event_pos[data['strand']=='-',5]
        data.loc[data['strand']=='-','exon_alt1_start']=event_pos[data['strand']=='-',0]+1
        data.loc[data['strand']=='-','exon_alt1_end']=event_pos[data['strand']=='-',1]
        data.loc[data['strand']=='-','exon_alt2_start']=event_pos[data['strand']=='-',2]+1
        data.loc[data['strand']=='-','exon_alt2_end']=event_pos[data['strand']=='-',3]
        iso1_str[data['strand']=='-']=data[data['strand']=='-'].apply(lambda x: str(x['exon_alt2_end'])+'j'+str(x['exon_const_start']),axis=1)
        iso2_str[data['strand']=='-']=data[data['strand']=='-'].apply(lambda x: str(x['exon_alt1_end'])+'j'+str(x['exon_const_start']),axis=1)

if event_type=='alt_5prime':
    data['exon_const_start']=event_pos[:,0]+1
    data['exon_const_end']=event_pos[:,1]
    data['exon_alt1_start']=event_pos[:,4]+1
    data['exon_alt1_end']=event_pos[:,5]
    data['exon_alt2_start']=event_pos[:,6]+1
    data['exon_alt2_end']=event_pos[:,7]
    iso1_str=data.apply(lambda x: str(x['exon_const_end'])+'j'+str(x['exon_alt2_start']),axis=1)
    iso2_str=data.apply(lambda x: str(x['exon_const_end'])+'j'+str(x['exon_alt1_start']),axis=1)
    if sum(data['strand']=='+')>0:
        data.loc[data['strand']=='+','exon_const_start']=event_pos[data['strand']=='+',4]+1
        data.loc[data['strand']=='+','exon_const_end']=event_pos[data['strand']=='+',5]
        data.loc[data['strand']=='+','exon_alt1_start']=event_pos[data['strand']=='+',0]+1
        data.loc[data['strand']=='+','exon_alt1_end']=event_pos[data['strand']=='+',1]
        data.loc[data['strand']=='+','exon_alt2_start']=event_pos[data['strand']=='+',2]+1
        data.loc[data['strand']=='+','exon_alt2_end']=event_pos[data['strand']=='+',3]
        iso1_str[data['strand']=='+']=data[data['strand']=='+'].apply(lambda x: str(x['exon_alt2_end'])+'j'+str(x['exon_const_start']),axis=1)
        iso2_str[data['strand']=='+']=data[data['strand']=='+'].apply(lambda x: str(x['exon_alt1_end'])+'j'+str(x['exon_const_start']),axis=1)

if event_type=='exon_skip':
    data['exon_pre_end']=event_pos[:,1]
    data['exon_start']=event_pos[:,2]+1
    data['exon_end']=event_pos[:,3]
    data['exon_aft_start']=event_pos[:,4]+1
    iso1_str=data.apply(lambda x: str(x['exon_pre_end'])+'j'+str(x['exon_start'])+'_'+str(x['exon_end'])+'j'+str(x['exon_aft_start']),axis=1)
    iso2_str=data.apply(lambda x: str(x['exon_pre_end'])+'j'+str(x['exon_aft_start']),axis=1)

if event_type=='mult_exon_skip':
    data['exon_pre_end']=event_pos[:,1]
    data['exon1_start']=event_pos[:,2]+1
    data['exon1_end']=event_pos[:,3]
    data['exon2_start']=event_pos[:,4]+1
    data['exon2_end']=event_pos[:,5]
    data['exon_aft_start']=event_pos[:,6]+1
    iso1_str=data.apply(lambda x: str(x['exon_pre_end'])+'j'+str(x['exon1_start'])+'_'+str(x['exon1_end'])+'j'+str(x['exon2_start'])+'_'+str(x['exon2_end'])+'j'+str(x['exon_aft_start']),axis=1)
    iso2_str=data.apply(lambda x: str(x['exon_pre_end'])+'j'+str(x['exon_aft_start']),axis=1)

if event_type=='intron_retention':
    data['exon1_end']=event_pos[:,1]
    data['exon2_start']=event_pos[:,2]+1
    iso1_str='NONE'
    iso2_str=data.apply(lambda x: str(x['exon1_end'])+'j'+str(x['exon2_start']),axis=1)

if event_type=='mutex_exons':
    tmp_reshape = np.transpose(event_pos, (0,2,1))
    event_coord=np.reshape(tmp_reshape,(len(gene_idx),8))
    data['exon_pre_end']=event_coord[:,1]
    data['exon1_start']=event_coord[:,2]+1
    data['exon1_end']=event_coord[:,3]
    data['exon2_start']=event_coord[:,4]+1
    data['exon2_end']=event_coord[:,5]
    data['exon_aft_start']=event_coord[:,6]+1
    if spladder_ver==1:
        switch_idx=data['exon1_end']-data['exon1_start']>data['exon2_end']-data['exon2_start']
        data.loc[switch_idx,'exon1_start']=event_coord[switch_idx,4]+1
        data.loc[switch_idx,'exon1_end']=event_coord[switch_idx,5]
        data.loc[switch_idx,'exon2_start']=event_coord[switch_idx,2]+1
        data.loc[switch_idx,'exon2_end']=event_coord[switch_idx,3]

    iso1_str=data.apply(lambda x: str(x['exon_pre_end'])+'j'+str(x['exon1_start'])+'_'+str(x['exon1_end'])+'j'+str(x['exon_aft_start']),axis=1)
    iso2_str=data.apply(lambda x: str(x['exon_pre_end'])+'j'+str(x['exon2_start'])+'_'+str(x['exon2_end'])+'j'+str(x['exon_aft_start']),axis=1)

data['event_jid']=data['contig'] + 's' + data['strand'] + ':g.' + iso1_str + '>' + iso2_str + "[spl" + event_type + "]"
if spladder_ver==1 and (event_type=='alt_3prime' or event_type=='alt_5prime'):
    iso1=pd.DataFrame(np.transpose(f['iso2'][idx1,start_pos:end_pos]),columns=list(map(lambda x: x + '_iso1',sample_names[idx1])))
    iso2=pd.DataFrame(np.transpose(f['iso1'][idx1,start_pos:end_pos]),columns=list(map(lambda x: x + '_iso2',sample_names[idx1])))
else:
    iso1=pd.DataFrame(np.transpose(f['iso1'][idx1,start_pos:end_pos]),columns=list(map(lambda x: x + '_iso1',sample_names[idx1])))
    iso2=pd.DataFrame(np.transpose(f['iso2'][idx1,start_pos:end_pos]),columns=list(map(lambda x: x + '_iso2',sample_names[idx1])))
if spladder_ver==1 and event_type=='mutex_exons':
    temp1=iso1.loc[switch_idx,:].copy()
    temp2=iso2.loc[switch_idx,:].copy()
    temp1.columns=iso2.columns
    temp2.columns=iso1.columns
    iso1.loc[switch_idx,:]=temp2
    iso2.loc[switch_idx,:]=temp1

data=pd.concat([data,iso1,iso2],axis=1)
data.to_csv(outname,index=False)
