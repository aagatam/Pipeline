import pyensembl
import Bio.SeqIO
import Bio.Seq
import pandas as pd
import sys
import re
import math
import os
import io
from Bio import pairwise2
import utils as bb
import warnings

def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()

print(sys.version_info)
print(pd._version)
print(bb.__file__)

events_file=sys.argv[1]
event_type=sys.argv[2]
aapad=int(sys.argv[3])
out_name=sys.argv[4]
ensemble_release=int(sys.argv[5])
ref_fasta=sys.argv[6]

ensembl=pyensembl.EnsemblRelease(ensemble_release,'mus_musculus')
ref=Bio.SeqIO.to_dict(Bio.SeqIO.parse(ref_fasta,'fasta'))

###read coordinates table
if event_type=="intron_retention":
 col_types={"contig":str,"event_id":str,"gene":str,"exon1_end":int,"exon2_start":int,"event_jid":str}
elif event_type=="exon_skip":
 col_types={"contig":str,"event_id":str,"gene":str,"exon_pre_end":int,"exon_start":int,"exon_end":int,"exon_aft_start":int,"event_jid":str}
elif event_type=="mult_exon_skip":
 col_types={"contig":str,"event_id":str,"gene":str,"exon_pre_end":int,"exon1_start":int,"exon1_end":int,"exon2_start":int,"exon2_end":int,"exon_aft_start":int,"event_jid":str}
elif event_type=="mutex_exons":
 col_types={"contig":str,"event_id":str,"gene":str,"exon_pre_end":int,"exon1_start":int,"exon1_end":int,"exon2_start":int,"exon2_end":int,"exon_aft_start":int,"event_jid":str}
elif event_type=="alt_3prime" or event_type=="alt_5prime":
 col_types={"contig":str,"event_id":str,"gene":str,"strand":str,"exon_const_start":int,"exon_const_end":int,"exon_alt1_start":int,"exon_alt1_end":int,"exon_alt2_start":int,"exon_alt2_end":int,"event_jid":str}
elif event_type=="ALL":
 col_types={"contig":str,"event_id":str,"gene":str,"strand":str,"event_jid":str}
else:
 print("ERROR: invalid event type " + event_type + ".\nValid event types are IR, ES, MULT, MUT, A3 or A5")
 sys.exit(1)

events_table=pd.read_csv(events_file,usecols=col_types.keys(),dtype=col_types)
events_table=events_table.assign(coding_transcript_effect=None,top_effect_transcript=None,effect_type=None,iso1_pc=None,iso2_pc=None,other_pc=None,iso1_nc=None,iso2_nc=None,other_nc=None)
topEffect_list=[]
wt_file='.'.join([out_name,"refSeq.fasta"])
if os.path.exists(wt_file):
 os.remove(wt_file)
wt_fasta=open(wt_file,"a+")
novel_file='.'.join([out_name,"altSeq.fasta"])
novel_file_change='.'.join([out_name,"changeSeq.fasta"])
novel_file_change_sub='.'.join([out_name,"changeSeq_subst.fasta"])
if os.path.exists(novel_file):
 os.remove(novel_file)
novel_fasta=open(novel_file,"a+")
novel_fasta_change=open(novel_file_change,"a+")
novel_fasta_change_sub=open(novel_file_change_sub,"a+")
peptides_file='.'.join([out_name,"peptides.csv"])
if os.path.exists(peptides_file):
 os.remove(peptides_file)
peptides_table=open(peptides_file,"a+")
effectCol=["altPept","event_id","effectId","gene","orf_effect","aa_effect","topEffect","refPept","refIsoform","sourceName","insSeqLen","refSeqLen","delSeqLen","refSeqPos","altSeqPos"]
effectDF=pd.DataFrame(columns=effectCol)
effectDF.to_csv(peptides_table,header=True,index=False)

## for each event
###find matching transcripts
###determine effect_type
for index,row in events_table.iterrows():
 effectDF=pd.DataFrame(columns=effectCol)
 wt_list=[]
 novel_list=[]
 novel_list_change = []
 novel_list_change_sub = []
 event_coords=bb.jid_to_coords(events_table.loc[index,"event_jid"])
 #event_coords=bb.get_event_coords(row,event_type)
 iso1_str=events_table.loc[index,"event_jid"].split("g.")[1].split(">")[0]
 iso2_str=events_table.loc[index,"event_jid"].split("g.")[1].split(">")[1].split("[")[0]
 transcript_table=bb.find_matching_transcripts(ensembl,row["gene"],event_coords)
 events_table.loc[index,"iso1_pc"]='|'.join(transcript_table.index.values[(transcript_table.coding==True) & (transcript_table.matching_isoform=="iso1")])
 events_table.loc[index,"iso2_pc"]='|'.join(transcript_table.index.values[(transcript_table.coding==True) & (transcript_table.matching_isoform=="iso2")])
 events_table.loc[index,"other_pc"]='|'.join(transcript_table.index.values[(transcript_table.coding==True) & (transcript_table.matching_isoform=="none")])
 events_table.loc[index,"iso1_nc"]='|'.join(transcript_table.index.values[(transcript_table.coding==False) & (transcript_table.matching_isoform=="iso1")])
 events_table.loc[index,"iso2_nc"]='|'.join(transcript_table.index.values[(transcript_table.coding==False) & (transcript_table.matching_isoform=="iso2")])
 events_table.loc[index,"other_nc"]='|'.join(transcript_table.index.values[(transcript_table.coding==False) & (transcript_table.matching_isoform=="none")])
 if events_table.loc[index,"iso1_pc"]!='' and events_table.loc[index,"iso2_pc"]!='':
  events_table.loc[index,"coding_transcript_effect"]="IsoformSwitch"
 elif events_table.loc[index,"iso1_pc"]!='':
  events_table.loc[index,"coding_transcript_effect"]="NovelIso2"
 elif events_table.loc[index,"iso2_pc"]!='':
  events_table.loc[index,"coding_transcript_effect"]="NovelIso1"
 elif events_table.loc[index,"other_pc"]!='':
  events_table.loc[index,"coding_transcript_effect"]="Unknown"
 else:
  events_table.loc[index,"coding_transcript_effect"]="NonCoding"
 #print(events_table.loc[index,])
 top_coding_effect=''
 top_effectId=''
 top_effect_type=''
 top_effect_score=-5
 max_wt_seq=-1
 for tid,transcript_info in transcript_table.loc[(transcript_table.coding==True) & ((transcript_table.matching_isoform=="iso1") | (transcript_table.matching_isoform=="iso2"))].iterrows():
  #print(tid)
  curr_effect_score=-5
  transcript=ensembl.transcript_by_id(tid)
  wt_seq=Bio.Seq.translate(transcript.coding_sequence)
  coding_coord=pd.DataFrame(transcript.coding_sequence_position_ranges,columns=['start','end'])
  coding_coord=coding_coord.append(pd.Series({'start':min(transcript.stop_codon_positions),'end':max(transcript.stop_codon_positions)}),ignore_index=True).sort_values(by="start")
  new_coord=bb.get_new_coord(transcript_info["matching_isoform"],event_coords,coding_coord)
   #print(novel_junc)
   #print(coding_coord)
   #print(overlap)
   #print(new_coord)
  (mut_seq,coding_effect)=bb.make_prot_from_coord(transcript,new_coord,ref)
   #print(mut_seq)
  (wt_diff_pos,mut_diff_pos)=bb.find_seq_diff(wt_seq,mut_seq)
   #print(wt_diff_pos)
   #print(mut_diff_pos)
  desc=' '.join(["GN=" + transcript.gene.gene_name,"ID=" + events_table.loc[index,"event_id"] + "_" + transcript_info["matching_isoform"],"JID=" + events_table.loc[index,"event_jid"],"PC=" + str(wt_diff_pos[0]) + "-" + str(wt_diff_pos[1])])
  wt_list.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id,description=desc,seq=Bio.Seq.Seq(wt_seq)))
  if transcript_info["matching_isoform"]=="iso1":
   desc=' '.join(["GN=" + transcript.gene.gene_name,"ID=" + events_table.loc[index,"event_id"] + "_iso2","JID=" + events_table.loc[index,"event_jid"],"PC=" + str(mut_diff_pos[0]) + "-" + str(mut_diff_pos[1])])
   #novel_list.append(Bio.SeqRecord.SeqRecord(id= "en|" +transcript.id + ":g." + iso1_str + ">" + iso2_str,description=desc,seq=Bio.Seq.Seq(mut_seq)))
  else:
   desc=' '.join(["GN=" + transcript.gene.gene_name,"ID=" + events_table.loc[index,"event_id"] + "_iso1","JID=" + events_table.loc[index,"event_jid"],"PC=" + str(mut_diff_pos[0]) + "-" + str(mut_diff_pos[1])])
   #novel_list.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso2_str + ">" + iso1_str,description=desc,seq=Bio.Seq.Seq(mut_seq)))
#"mutPept","event_id","event_jid","gene_name","transcript_id","orf_effect","topEffect","wtPept","sourceName","novelSeqLen","wtSeqLen","delSeqLen"]
  if not math.isnan(mut_diff_pos[0]):
   mut_pept=mut_seq[max(mut_diff_pos[0]-aapad,0):mut_diff_pos[1]+aapad]
   wt_pept=wt_seq[max(mut_diff_pos[0]-aapad,0):wt_diff_pos[1]+aapad]
   novel_seq_len=mut_diff_pos[1]-mut_diff_pos[0]
   del_seq_len=wt_diff_pos[1]-wt_diff_pos[0]
  else:
   mut_pept=''
   wt_pept=''
   novel_seq_len=float('nan')
   del_seq_len=float('nan')
  if (coding_effect=="InFrame") and math.isnan(mut_diff_pos [0]):
   aa_effect="Silent"
   curr_effect_score=-4
  elif del_seq_len==len(wt_seq)-1:
   aa_effect="ProteinLoss"
   curr_effect_score=-2
   if transcript_info["matching_isoform"] == "iso1":
    desc = ' '.join(["GN=" + transcript.gene.gene_name, "ID=" + events_table.loc[index, "event_id"] + "_iso2", "Eff=" + aa_effect,
                     "JID=" + events_table.loc[index, "event_jid"],
                     "PC=" + str(mut_diff_pos[0]) + "-" + str(mut_diff_pos[1])])
    novel_list_change.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso1_str + ">" + iso2_str, description=desc,
                             seq=Bio.Seq.Seq(wt_pept)))
    desc = ' '.join(["GN=" + transcript.gene.gene_name, "ID=" + events_table.loc[index, "event_id"] + "_iso2", "Eff=" + aa_effect,
                     "JID=" + events_table.loc[index, "event_jid"],
                     "PC=" + str(mut_diff_pos[0]) + "-" + str(mut_diff_pos[1])])
    novel_list.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso1_str + ">" + iso2_str, description=desc,
                             seq=Bio.Seq.Seq(wt_seq)))
   else:
    desc = ' '.join(["GN=" + transcript.gene.gene_name, "ID=" + events_table.loc[index, "event_id"] + "_iso1", "Eff=" + aa_effect,
                     "JID=" + events_table.loc[index, "event_jid"],
                     "PC=" + str(mut_diff_pos[0]) + "-" + str(mut_diff_pos[1])])

    novel_list_change.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso2_str + ">" + iso1_str, description=desc,
                             seq=Bio.Seq.Seq(wt_pept)))
    desc = ' '.join(["GN=" + transcript.gene.gene_name, "ID=" + events_table.loc[index, "event_id"] + "_iso1", "Eff=" + aa_effect,
                     "JID=" + events_table.loc[index, "event_jid"],
                     "PC=" + str(mut_diff_pos[0]) + "-" + str(mut_diff_pos[1])])

    novel_list.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso2_str + ">" + iso1_str, description=desc,
                             seq=Bio.Seq.Seq(wt_seq)))

  elif (novel_seq_len>0) and (del_seq_len>0):
   aa_effect="Substitution"
   curr_effect_score=novel_seq_len
   if len(mut_pept)!=0:
    if transcript_info["matching_isoform"] == "iso1":
     desc = ' '.join(
      ["GN=" + transcript.gene.gene_name, "ID=" + events_table.loc[index, "event_id"] + "_iso2", "Eff=" + aa_effect,
       "JID=" + events_table.loc[index, "event_jid"],
       "PC=" + str(mut_diff_pos[0]) + "-" + str(mut_diff_pos[1])])
     novel_list_change.append(
      Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso1_str + ">" + iso2_str, description=desc,
                              seq=Bio.Seq.Seq(mut_pept)))
     novel_list_change_sub.append(
      Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso1_str + ">" + iso2_str, description=desc,
                              seq=Bio.Seq.Seq(wt_pept)))
     novel_list.append(
      Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso1_str + ">" + iso2_str, description=desc,
                              seq=Bio.Seq.Seq(mut_seq)))
    else:
     desc = ' '.join(
      ["GN=" + transcript.gene.gene_name, "ID=" + events_table.loc[index, "event_id"] + "_iso1", "Eff=" + aa_effect,
       "JID=" + events_table.loc[index, "event_jid"],
       "PC=" + str(mut_diff_pos[0]) + "-" + str(mut_diff_pos[1])])

     novel_list_change.append(
      Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso2_str + ">" + iso1_str, description=desc,
                              seq=Bio.Seq.Seq(mut_pept)))
     novel_list_change_sub.append(
      Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso2_str + ">" + iso1_str, description=desc,
                              seq=Bio.Seq.Seq(wt_pept)))
     novel_list.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso2_str + ">" + iso1_str, description=desc,
                              seq=Bio.Seq.Seq(mut_seq)))
  elif (novel_seq_len>0):
   aa_effect="Insertion"
   curr_effect_score=novel_seq_len
   if len(mut_pept)!=0:
    if transcript_info["matching_isoform"] == "iso1":
     desc = ' '.join(
      ["GN=" + transcript.gene.gene_name, "ID=" + events_table.loc[index, "event_id"] + "_iso2", "Eff=" + aa_effect,
       "JID=" + events_table.loc[index, "event_jid"],
       "PC=" + str(mut_diff_pos[0]) + "-" + str(mut_diff_pos[1])])
     novel_list_change.append(
      Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso1_str + ">" + iso2_str, description=desc,
                              seq=Bio.Seq.Seq(mut_pept)))
     novel_list.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso1_str + ">" + iso2_str, description=desc,
                              seq=Bio.Seq.Seq(mut_seq)))
    else:
     desc = ' '.join(
      ["GN=" + transcript.gene.gene_name, "ID=" + events_table.loc[index, "event_id"] + "_iso1", "Eff=" + aa_effect,
       "JID=" + events_table.loc[index, "event_jid"],
       "PC=" + str(mut_diff_pos[0]) + "-" + str(mut_diff_pos[1])])

     novel_list_change.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso2_str + ">" + iso1_str, description=desc,
                              seq=Bio.Seq.Seq(mut_pept)))
     novel_list.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso2_str + ">" + iso1_str, description=desc,
                              seq=Bio.Seq.Seq(mut_seq)))
  elif (del_seq_len>0) and (wt_diff_pos[1]==len(wt_seq)-1):
   aa_effect="Truncation"
   curr_effect_score=-3
   if transcript_info["matching_isoform"] == "iso1":
    desc = ' '.join(["GN=" + transcript.gene.gene_name, "ID=" + events_table.loc[index, "event_id"] + "_iso2", "Eff=" + aa_effect,
                     "JID=" + events_table.loc[index, "event_jid"],
                     "PC=" + str(mut_diff_pos[0]) + "-" + str(mut_diff_pos[1])])
    novel_list_change.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso1_str + ">" + iso2_str, description=desc,
                             seq=Bio.Seq.Seq(wt_pept)))
    novel_list.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso1_str + ">" + iso2_str, description=desc,
                             seq=Bio.Seq.Seq(mut_seq)))
   else:
    desc = ' '.join(["GN=" + transcript.gene.gene_name, "ID=" + events_table.loc[index, "event_id"] + "_iso1", "Eff=" + aa_effect,
                     "JID=" + events_table.loc[index, "event_jid"],
                     "PC=" + str(mut_diff_pos[0]) + "-" + str(mut_diff_pos[1])])

    novel_list_change.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso2_str + ">" + iso1_str, description=desc,
                             seq=Bio.Seq.Seq(wt_pept)))
    novel_list.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso2_str + ">" + iso1_str, description=desc,
                             seq=Bio.Seq.Seq(mut_seq)))
  elif (del_seq_len>0):
   aa_effect="Deletion"
   curr_effect_score=-1
   if transcript_info["matching_isoform"] == "iso1":
    desc = ' '.join(["GN=" + transcript.gene.gene_name, "ID=" + events_table.loc[index, "event_id"] + "_iso2", "Eff=" + aa_effect,
                     "JID=" + events_table.loc[index, "event_jid"],
                     "PC=" + str(mut_diff_pos[0]) + "-" + str(mut_diff_pos[1])])
    novel_list_change.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso1_str + ">" + iso2_str, description=desc,
                             seq=Bio.Seq.Seq(wt_pept)))
    novel_list.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso1_str + ">" + iso2_str, description=desc,
                             seq=Bio.Seq.Seq(mut_seq)))
   else:
    desc = ' '.join(["GN=" + transcript.gene.gene_name, "ID=" + events_table.loc[index, "event_id"] + "_iso1", "Eff=" + aa_effect,
                     "JID=" + events_table.loc[index, "event_jid"],
                     "PC=" + str(mut_diff_pos[0]) + "-" + str(mut_diff_pos[1])])

    novel_list_change.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso2_str + ">" + iso1_str, description=desc,
                             seq=Bio.Seq.Seq(wt_pept)))
    novel_list.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso2_str + ">" + iso1_str, description=desc,
                             seq=Bio.Seq.Seq(mut_seq)))
  else:
   aa_effect="Unknown"
  effectId=tid + "_" +  events_table.loc[index,"event_jid"]
  if (curr_effect_score>top_effect_score) or ((curr_effect_score==top_effect_score) and (max_wt_seq>len(wt_seq))):
   top_coding_effect=aa_effect + "_" + coding_effect
   top_effectId=effectId
   top_effect_score=curr_effect_score
   max_wt_seq=len(wt_seq)
   curr_effect=pd.Series({"altPept": mut_pept, "event_id": events_table.loc[index,"event_id"], "effectId": effectId, "gene": transcript.gene.gene_name,"orf_effect": coding_effect,"aa_effect": aa_effect,"refPept": wt_pept,"refIsoform": transcript_info["matching_isoform"],"insSeqLen": novel_seq_len,"refSeqLen": len(wt_seq)-1,"delSeqLen":del_seq_len,"refSeqPos":str(wt_diff_pos[0]+1) + '-' + str(wt_diff_pos[1]+1),"altSeqPos":str(mut_diff_pos[0]+1) + '-' + str(mut_diff_pos[1]+1)})
   #print(curr_effect)
   effectDF=effectDF.append(curr_effect,ignore_index=True)

 out_handle=io.StringIO()
 Bio.SeqIO.write(wt_list,out_handle,"fasta")
 wt_fasta.write(out_handle.getvalue())
 out_handle=io.StringIO()
 Bio.SeqIO.write(novel_list,out_handle,"fasta")
 novel_fasta.write(out_handle.getvalue())
 out_handle = io.StringIO()
 Bio.SeqIO.write(novel_list_change, out_handle, "fasta")
 novel_fasta_change.write(out_handle.getvalue())
 out_handle = io.StringIO()
 Bio.SeqIO.write(novel_list_change_sub, out_handle, "fasta")
 novel_fasta_change_sub.write(out_handle.getvalue())
 effectDF.topEffect=effectDF.effectId==top_effectId
 effectDF.to_csv(peptides_table,header=False,index=False)
 events_table.loc[index,"effect_type"]=top_coding_effect
 events_table.loc[index,"top_effect_transcript"]=top_effectId.split('_')[0]


events_table.to_csv('.'.join([out_name,'effects.csv']))
wt_fasta.close()
novel_fasta.close()
novel_fasta_change.close()
peptides_table.close()
#Bio.SeqIO.write(wt_list,'.'.join([out_name,event_type,"wtSeq.fasta"]),"fasta")
#Bio.SeqIO.write(novel_list,'.'.join([out_name,event_type,"novelSeq.fasta"]),"fasta")
#### for each isoform switch
#####add seq to fastas

#### for each iso1/iso2 seq
#####make new coord
#####make new seq
#####make effect id
#####compare seq
#####classify effect
#####add seq to fastas
#####get mut/wt peptides
#####score effect

### find top effect
