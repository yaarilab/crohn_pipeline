$HOSTNAME = ""
params.outdir = 'results'  

// Add for each process an option to change the parameters. Default is the set params
//* params.edit_IgBlast_params =  "no"  //* @dropdown @options:"yes","no"  @show_settings:"IgBlast"
//* params.edit_MakeDb_igblast_params =  "no"  //* @dropdown @options:"yes","no" @show_settings:"MakeDb_igblast"
//* params.edit_CreateGermlines_params =  "no"  //* @dropdown @options:"yes","no" @show_settings:"CreateGermlines"
//* params.edit_DefineClones_params =  "no"  //* @dropdown @options:"yes","no" @show_settings:"DefineClones"
//* params.edit_change_names_germ_pass_params =  "no"  //* @dropdown @options:"yes","no" @show_settings:"change_names_germ_pass"
//* params.edit_change_names_clone_pass_params =  "no"  //* @dropdown @options:"yes","no" @show_settings:"change_names_clone_pass"
//* params.edit_crohn_plot_params =  "no"  //* @dropdown @options:"yes","no" @show_settings:"crohn_analysis"
//* params.edit_pipeAIRR_figure2_panelC_params =  "no"  //* @dropdown @options:"yes","no" @show_settings:"pipeAIRR_figure2_panelC"

//* autofill
if ($HOSTNAME == "default"){
    $DOCKER_IMAGE = "immcantation/suite:4.3.0"
    $DOCKER_OPTIONS = "-v /work:/work"

}

//* platform
if ($HOSTNAME == "ig03.lnx.biu.ac.il"){
    $DOCKER_IMAGE = "immcantation/suite:4.3.0"
    $DOCKER_OPTIONS = "-v /work:/work"
	$CPU  = 48
    $MEMORY = 300 
}
//* platform


//* autofill

run_dir = params.outdir.replace("report", "run")

if((params.edit_IgBlast_params && (params.edit_IgBlast_params == "no"))){
    // Process Parameters for params.IgBlast:
    params.IgBlast.num_threads = "60" 
	params.IgBlast.db_v = "${run_dir}/igblast_reference/human_gl_IGHV_F+ORF+in-frame_P_gaps_new.fasta" 
	params.IgBlast.db_d = "/${run_dir}/igblast_reference/human_gl_IGHD_F+ORF+in-frame_P.fasta" 
	params.IgBlast.db_j = "${run_dir}/igblast_reference/human_gl_IGHJ_F+ORF+in-frame_P.fasta"
	params.IgBlast.num_alignments_V = "10"
	params.IgBlast.domain_system = "imgt"
	params.IgBlast.auxiliary_data = "/usr/local/share/igblast/optional_file/human_gl.aux"
	params.IgBlast.outfmt = "MakeDb"
}

if((params.edit_MakeDb_igblast_params && (params.edit_MakeDb_igblast_params == "no"))){
    // Process Parameters for params.MakeDb_igblast:
    params.MakeDb_igblast.failed = "false"
	params.MakeDb_igblast.v_file = "${run_dir}/igblast_reference/human_gl_IGHV_F+ORF+in-frame_P_w_gaps_new.fasta" 
	params.MakeDb_igblast.d_file = "${run_dir}/igblast_reference/human_gl_IGHD_F+ORF+in-frame_P.fasta" 
	params.MakeDb_igblast.j_file = "${run_dir}/igblast_reference/human_gl_IGHJ_F+ORF+in-frame_P.fasta"
    params.MakeDb_igblast.extended = "true"
    params.MakeDb_igblast.format  =  "airr" 
    params.MakeDb_igblast.regions = "default"
	params.MakeDb_igblast.asisid = "false"
	params.MakeDb_igblast.asiscalls = "false"
	params.MakeDb_igblast.inferjunction = "fasle"
	params.MakeDb_igblast.partial = "false"
}

if((params.edit_CreateGermlines_params && (params.edit_CreateGermlines_params == "no"))){
	// Process Parameters for params.CreateGermlines:
	params.CreateGermlines.failed = "false" 
	params.CreateGermlines.format =  "airr"
	params.CreateGermlines.references = "${run_dir}/germlines_reference/imgt_human_IGHV.fasta  ${run_dir}/germlines_reference/imgt_human_IGHD.fasta   ${run_dir}/germlines_reference/imgt_human_IGHJ.fasta"
	params.CreateGermlines.g = "dmask" 
	params.CreateGermlines.cloned = "false"
	params.CreateGermlines.seq_field = "" 
	params.CreateGermlines.v_field = "" 
	params.CreateGermlines.d_field = "" 
	params.CreateGermlines.j_field = "" 
	params.CreateGermlines.clone_field = ""
}

if((params.edit_DefineClones_params && (params.edit_DefineClones_params == "no"))){
	// Process Parameters for params.DefineClones:
	params.DefineClones.failed = "false" 
	params.DefineClones.format =  "airr"
	params.DefineClones.seq_field = ""  
	params.DefineClones.v_field = ""  
	params.DefineClones.d_field = "" 
	params.DefineClones.j_field = ""  
	params.DefineClones.group_fields = "" 
	
	params.DefineClones.mode = "gene" 
	params.DefineClones.dist = "0.2"  
	params.DefineClones.norm = "len"
	params.DefineClones.act = "set" 
	params.DefineClones.model = "ham" 
	params.DefineClones.sym = "avg"
	params.DefineClones.link = "single" 
	params.DefineClones.maxmiss = "0"
}


if((params.edit_change_names_germ_pass_params && (params.edit_change_names_germ_pass_params == "no"))){
	 // Process Parameters for params.change_names_germ_pass:
	params.change_names_germ_pass.tsv_file = "${run_dir}/tables/name_table.tsv"
}

if((params.edit_change_names_clone_pass_params && (params.edit_change_names_clone_pass_params == "no"))){
	 // Process Parameters for params.change_names_clone_pass:
	params.change_names_clone_pass.tsv_file = "${run_dir}/tables/name_table.tsv"
}

if((params.edit_crohn_plot_params && (params.edit_crohn_plot_params == "no"))){
	 // Process Parameters for params.crohn_analysis:
	params.crohn_analysis.sup1_table = "${run_dir}/sup1_table.tsv"
}
if((params.edit_pipeAIRR_figure2_panelC_params && (params.edit_pipeAIRR_figure2_panelC_params == "no"))){
	 // Process Parameters for params.pipeAIRR_figure2_panelC:
	params.pipeAIRR_figure2_panelC.original_v_usage = "${run_dir}/v_usage_safra_et_al.tsv"
}

if (!params.inputparam){params.inputparam = ""} 
if (!params.mate){params.mate = ""} 

if (params.inputparam){
Channel
	.fromFilePairs( params.inputparam , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.inputparam}" }
	.set{g_9_reads_g_13}
 } else {  
	g_9_reads_g_13 = Channel.empty()
 }

Channel.value(params.mate).set{g_14_mate_g_13}


process fastqTofasta {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.fasta$/) "outputparam/$filename"}
input:
 set val(name),file(reads) from g_9_reads_g_13
 val mate from g_14_mate_g_13

output:
 set val(name),file("*.fasta")  into g_13_fastaFile0_g_10

script:
	
readArray = reads.toString().split(' ')	

if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	R1n = R1.replace('.fastq','')
	R2n = R2.replace('.fastq','')
	
	"""
	 awk 'BEGIN {n_seq = 0;count = 0;} {header = ">seq" ++n_seq;} {if (n_seq % 4 == 1) {l=">" ++count ;print l; ;} else if (n_seq % 4 == 2) {print;}}' ${R1n}.fastq > ${R1n}.fasta
	 awk 'BEGIN {n_seq = 0;count = 0;} {header = ">seq" ++n_seq;} {if (n_seq % 4 == 1) {l=">" ++count ;print l; ;} else if (n_seq % 4 == 2) {print;}}' ${R2n}.fastq > ${R2n}.fasta
	"""
	
}else{

	"""
	 awk 'BEGIN {n_seq = 0;count = 0;} {header = ">seq" ++n_seq;} {if (n_seq % 4 == 1) {l=">" ++count ;print l; ;} else if (n_seq % 4 == 2) {print;}}' ${name}.fastq > ${name}.fasta
	"""
}
}


process IgBlast {

input:
 set val(name),file(fastaFile) from g_13_fastaFile0_g_10

output:
 file "*igblast*"  into g_10_logFile0_g_28
 set val(name),file(fastaFile)  into g_10_fastaFile1_g_28

script:
num_threads = params.IgBlast.num_threads
db_v = params.IgBlast.db_v
db_d = params.IgBlast.db_d
db_j = params.IgBlast.db_j
num_alignments_V = params.IgBlast.num_alignments_V
domain_system = params.IgBlast.domain_system
auxiliary_data = params.IgBlast.auxiliary_data
outfmt = params.IgBlast.outfmt

outfile = (outfmt=="MakeDb") ? name+"_igblast.out" : name+"_igblast.tsv"
outfmt = (outfmt=="MakeDb") ? "'7 std qseq sseq btop'" : "19"

"""
mkdir db_blast

export IGDATA=/usr/local/share/igblast

makeblastdb -parse_seqids -dbtype nucl -in ${db_v} -out db_blast/V

makeblastdb -parse_seqids -dbtype nucl -in ${db_d} -out db_blast/D

makeblastdb -parse_seqids -dbtype nucl -in ${db_j} -out db_blast/J

igblastn -query ${fastaFile} -germline_db_V db_blast/V -germline_db_D db_blast/D  -germline_db_J db_blast/J -num_alignments_V ${num_alignments_V} -domain_system ${domain_system} -auxiliary_data ${auxiliary_data}  -outfmt ${outfmt}  -num_threads ${num_threads} -out ${outfile}>>  ${name}_igblast.out  
"""

}


process MakeDb_igblast {

input:
 set val(name),file(fastaFile) from g_10_fastaFile1_g_28
 file igblastfile from g_10_logFile0_g_28

output:
 set val(name),file("*_db-pass.tsv")  into g_28_outputFileTSV0_g_17

script:
failed = params.MakeDb_igblast.failed
v_file = params.MakeDb_igblast.v_file
d_file = params.MakeDb_igblast.d_file
j_file = params.MakeDb_igblast.j_file
format = params.MakeDb_igblast.format
regions = params.MakeDb_igblast.regions
extended = params.MakeDb_igblast.extended
asisid = params.MakeDb_igblast.asisid
asiscalls = params.MakeDb_igblast.asiscalls
inferjunction = params.MakeDb_igblast.inferjunction
partial = params.MakeDb_igblast.partial


readArray = fastaFile.toString().split(' ')
file = readArray.grep(~/.*.fasta*/)[0]
failed = (failed=="true") ? "--failed" : ""

format = (format=="changeo") ? "--format changeo" : ""
extended = (extended=="true") ? "--extended" : ""


regions = (regions=="rhesus-igl") ? "--regions rhesus-igl" : ""
asisid = (asisid=="true") ? "--asis-id" : ""
asiscalls = (asiscalls=="true") ? "--asis-calls" : ""
inferjunction = (inferjunction=="true") ? "--infer-junction" : ""
partial = (partial=="true") ? "--partial" : ""

"""
MakeDb.py igblast -s ${file} -i ${igblastfile} -r ${v_file} ${d_file} ${j_file} ${format} ${extended} ${regions} ${asisid} ${asiscalls} ${inferjunction} ${partial} --log MD_${name}.log ${failed}
"""


}


process CreateGermlines {

input:
 set val(name),file(reads) from g_28_outputFileTSV0_g_17

output:
 set val(name),file("*_germ-pass.tsv")  into g_17_outputFileTSV0_g_18, g_17_outputFileTSV0_g_26

script:
failed = params.CreateGermlines.failed
format = params.CreateGermlines.format

references = params.CreateGermlines.references
g = params.CreateGermlines.g
cloned = params.CreateGermlines.cloned
seq_field = params.CreateGermlines.seq_field
v_field = params.CreateGermlines.v_field
d_field = params.CreateGermlines.d_field
j_field = params.CreateGermlines.j_field
clone_field = params.CreateGermlines.clone_field


failed = (failed=="true") ? "--failed" : ""
format = (format=="airr") ? "": "--format changeo"
g = (g=="dmask") ? "" : "-g ${g}"
cloned = (cloned=="false") ? "" : "--cloned"


	
v_field = (v_field=="") ? "" : "--vf ${v_field}"
d_field = (d_field=="") ? "" : "--df ${d_field}"
j_field = (j_field=="") ? "" : "--jf ${j_field}"
seq_field = (seq_field=="") ? "" : "--sf ${seq_field}"



"""
echo ${format}

CreateGermlines.py -d ${reads} -r ${references} ${failed} ${format} ${g} ${cloned} ${v_field} ${d_field} ${j_field} ${seq_field} ${clone_field} --log CG_${name}.log 
"""



}


process DefineClones {

input:
 set val(name),file(reads) from g_17_outputFileTSV0_g_18

output:
 set val(name),file("*_clone-pass.tsv")  into g_18_outputFileTSV0_g_31

script:
failed = params.DefineClones.failed
format = params.DefineClones.format
seq_field = params.DefineClones.seq_field
v_field = params.DefineClones.v_field
d_field = params.DefineClones.d_field
j_field = params.DefineClones.j_field
group_fields = params.DefineClones.group_fields

mode = params.DefineClones.mode
dist = params.DefineClones.dist
norm = params.DefineClones.norm
act = params.DefineClones.act
model = params.DefineClones.model
sym = params.DefineClones.sym
link = params.DefineClones.link
maxmiss = params.DefineClones.maxmiss

failed = (failed=="true") ? "--failed" : ""
format = (format=="airr") ? "--format airr": "--format changeo"
group_fields = (group_fields=="") ? "" : "--gf ${group_fields}"
v_field = (v_field=="") ? "" : "--vf ${v_field}"
d_field = (d_field=="") ? "" : "--df ${d_field}"
j_field = (j_field=="") ? "" : "--jf ${j_field}"
seq_field = (seq_field=="") ? "" : "--sf ${seq_field}"


mode = (mode=="gene") ? "" : "--mode ${mode}"
norm = (norm=="len") ? "" : "--norn ${norm}"
act = (act=="set") ? "" : "--act ${act}"
model = (model=="ham") ? "" : "--model {model}"
sym = (sym=="avg") ? "" : "--sym {sym}"
link = (link=="single") ? "" : " --link ${link}"
    
	
"""
DefineClones.py -d ${reads} ${failed} ${format} ${v_field} ${d_field} ${j_field} ${seq_field} ${group_fields} ${mode} ${act} ${model} --dist ${dist} ${norm} ${sym} ${link} --maxmiss ${maxmiss} --log DF_.log  
"""



}


process change_names_clone_pass {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_clone-pass.tsv$/) "files/$filename"}
input:
 set val(name),file(clonePass) from g_18_outputFileTSV0_g_31

output:
 file "*_clone-pass.tsv"  into g_31_outputFileTSV0_g_24

script:
tsv_file = params.change_names_clone_pass.tsv_file


readArrayC = clonePass.toString().split(' ')	


"""
#!/usr/bin/env Rscript 
print("${tsv_file}")
d<-read.csv("${tsv_file}" ,sep="\t")
new_name = d[which(d[,1]=="${name}",),2]
file.rename("${readArrayC[0]}",paste(new_name ,"_germ-pass_clone-pass.tsv",sep=""))
print(new_name)
"""

}


process change_names_germ_pass {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_germ-pass.tsv$/) "files/$filename"}
input:
 set val(name),file(germPass) from g_17_outputFileTSV0_g_26

output:
 file "*_germ-pass.tsv"  into g_26_outputFileTSV0_g_24



script:
tsv_file = params.change_names_germ_pass.tsv_file

readArrayG = germPass.toString().split(' ')	
"""
#!/usr/bin/env Rscript 
print("${tsv_file}")
d<-read.csv("${tsv_file}" ,sep="\t")
new_name = d[which(d[,1]=="${name}",),2]
file.rename("${readArrayG[0]}",paste(new_name ,"_germ-pass.tsv",sep=""))
print(new_name)
"""



}


process crohn_analysis {

input:
 file gremPass from g_26_outputFileTSV0_g_24.collect()
 file clonePass from g_31_outputFileTSV0_g_24.collect()

output:
 file "*basics.csv"  into g_24_csvFile0_g_32
 file "*q.csv"  into g_24_csvFile1_g_32
 file "*V.csv"  into g_24_csvFile2_g_32, g_24_csvFile2_g_33
 file "*VJ.csv"  into g_24_csvFile3_g_32

script:

sup1_table = params.crohn_analysis.sup1_table


readArrayG = gremPass.toString().split(' ').toString()
readArrayG = readArrayG.replace('[','')
readArrayG = readArrayG.replace(']','')

readArrayC = clonePass.toString().split(' ').toString()
readArrayC = readArrayC.replace('[','')
readArrayC = readArrayC.replace(']','')


"""
#!/usr/bin/env Rscript 


dir="${params.outdir}"

xG<-strsplit("${readArrayG}",", ")
xC<-strsplit("${readArrayC}",", ")
fG=unlist(xG)
fC=unlist(xC)

fG <- paste0(paste(dir,"/files/", sep=""), fG)
fC <- paste0(paste(dir,"/files/", sep=""), fC)

f<-append(fG,fC)


library(data.table)
library(parallel)

p=mclapply(f,function(i){
  print(i)
  path<-paste(dir,"/files/",i, sep="")
  r=data.frame(id=gsub('_germ-pass.tsv','',gsub(paste(dir,"/files/", sep=""),'',i)))
  rep=read.delim(i)
  l=sort(nchar(rep[!is.na(rep[,"junction_length"]),][,"junction_aa"]))
  r=merge(r,data.frame(jl10=l[0.1*length(l)],
                       jl50=median(l),jl90=l[0.9*length(l)]))
  l=sort(rep[!is.na(rep[,"v_identity"]),][,"v_identity"])
  r=merge(r,data.frame(vi10=l[0.1*length(l)],
                       vi50=median(l),vi90=l[0.9*length(l)]))
                       
  
},mc.cores=50)

library(data.table)
z=data.frame(rbindlist(p))
db=read.csv("${sup1_table}", sep="\t")
db=db[db[,"Age"]<=18 & db[,"Clinical.state"]!='case ',]
db=db[!duplicated(db[,"Volunteer"]),]
x<-db[,c(1,3)]
colnames(x) = c("id", "stage")
z=merge(z,x,by='id')
z[,"stage"]=ifelse(grepl('case',z[,"stage"]),'case','control')
z["id"]=NULL
z[,4]=(1-z[,4]);z[,5]=(1-z[,5]);z[,6]=(1-z[,6])
basics=z
library(alakazam)
library(plyr)
p=mclapply(f,function(i){

  id=id=gsub('_germ-pass.tsv','',gsub(paste(dir,"/files/", sep=""),'',i))

  rep=read.delim(i)
  rep[,"v"]=getGene(rep[,"v_call"])
  rep[,"j"]=getGene(rep[,"j_call"])
  rep[,"VF"]=getFamily(rep[,"v_call"])
  rep[,"vj"]=paste(getGene(rep[,"v_call"]),getGene(rep[,"j_call"]))
  r=data.frame(table(rep[,"vj"]))
  r[,"Freq"]=r[,"Freq"]/sum(r[,"Freq"])
  r[,"id"]=id
  return(r)
},mc.cores=50)
p=data.frame(rbindlist(p))
m=merge(levels(as.factor(p[,"id"])),levels(as.factor(p[,"Var1"])))
m[,"p"]=paste(m[,"x"],m[,"y"])
colnames(m)=c('id','Var1','p')
p[,"p"]=paste(p[,"id"],p[,"Var1"])
m=merge(m,p,by='p',all=T)
m[,"p"]=NULL;m[,"Var1.y"]=NULL;m[,"id.y"]=NULL
colnames(m)=c('id','vj','freq')
m[,"freq"]=ifelse(is.na(m[,"freq"]),0,m[,"freq"])
z=merge(m,x,by='id')
z[,"stage"]=ifelse(grepl('case',z[,"stage"]),'case','control')

l=count(z,'vj','freq')
l=l[order(l[,"freq"],decreasing = T),]
l=l[1:50,]
z=z[z[,"vj"]%in%l[,"vj"],]
p=list()
j=1
for(i in unique(z[,"vj"])){
  p[[j*2]]=z[z[,"stage"]=='case'&z[,"vj"]==i,][,"freq"]
  p[[j*2-1]]=z[z[,"stage"]!='case'&z[,"vj"]==i,][,"freq"]
  print(t.test(p[[j*2]],p[[j*2-1]])["p.value"])
  j=j+1
}
names(p)=rep(unique(z[,"vj"]),each=2)
vj=p

p=mclapply(f,function(i){
  id=id=gsub('_germ-pass.tsv','',gsub(paste(dir,"/files/", sep=""),'',i))
  rep=read.delim(i)
  rep[,"v"]=getGene(rep[,"v_call"])
  rep[,"j"]=getGene(rep["j_call"])
  rep[,"VF"]=getFamily(rep[,"v_call"])
  rep[,"vj"]=paste(getGene(rep[,"v_call"]),getGene(rep[,"j_call"]))
  r=data.frame(table(rep[,"v"]))
  r[,"Freq"]=r[,"Freq"]/sum(r[,"Freq"])
  r[,"id"]=id
  return(r)
},mc.cores=50)
p=data.frame(rbindlist(p))
m=merge(levels(as.factor(p[,"id"])),levels(as.factor(p[,"Var1"])))
m[,"p"]=paste(m[,"x"],m[,"y"])
colnames(m)=c('id','Var1','p')
p[,"p"]=paste(p[,"id"],p[,"Var1"])
m=merge(m,p,by='p',all=T)
m[,"p"]=NULL;m[,"Var1.y"]=NULL;m[,"id.y"]=NULL
colnames(m)=c('id','vj','freq')
m[,"freq"]=ifelse(is.na(m[,"freq"]),0,m[,"freq"])
z=merge(m,x,by='id')
z[,"stage"]=ifelse(grepl('case',z[,"stage"]),'case','control')

l=count(z,'vj','freq')
l=l[order(l[,"freq"],decreasing = T),]
l=l[1:50,]
z=z[z[,"vj"]%in%l[,"vj"],]
p=list()
j=1
for(i in unique(z[,"vj"])){
  p[[j*2]]=z[z[,"stage"]=='case'&z[,"vj"]==i,][,"freq"]
  p[[j*2-1]]=z[z[,"stage"]!='case'&z[,"vj"]==i,][,"freq"]
  print(t.test(p[[j*2]],p[[j*2-1]])["p.value"])
  j=j+1
}
names(p)=rep(unique(z[,"vj"]),each=2)
v=p



library(shazam)
library(alakazam)
library(parallel)

p=mclapply(fC,function(i){
  a=read.delim(i)
  print(paste(i,nrow(a)))
  if(nrow(a>8000)){
    aa=data.frame(clone_id=a[!is.na(a[,"clone_id"]),][,"clone_id"])
    e=alphaDiversity(a,max_n=8000)
    e@diversity[,"group"]=gsub(paste(dir,"/files/", sep=""),'',
                           gsub('_germ-pass_clone-pass.tsv','',i))
    return(e@diversity)
  }else return(NULL)
},mc.cores=50)
library(data.table)
p2=rbindlist(p)

db=read.csv("${sup1_table}", sep="\t")
db=db[db[,"Age"]<=18 & db[,"Clinical.state"]!='case ',]
db=db[!duplicated(db[,"Volunteer"]),]
x<-db[,c(1,3)]
colnames(x) = c("id", "stage")
z=merge(z,x,by='id')
p=merge(p2,x,by.x='group',by.y='id')
p=data.frame(p)
p[,"stage"]=ifelse(grepl('case',p[,"stage"]),'case','control')
q=p[p[,"q"]%in%c(0,1,2,3,4),]

write.csv(basics,"basics.csv")
write.csv(q ,"q.csv")
write.csv(v,"V.csv")
write.csv(vj,"VJ.csv")

"""

}


process pipeAIRR_figure2_panelC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*pdf$/) "pdf/$filename"}
input:
 file v_usage from g_24_csvFile2_g_33

output:
 file "*pdf"  into g_33_outputFilePdf00

script:

original_v_usage = params.pipeAIRR_figure2_panelC.original_v_usage

"""
#!/usr/bin/env Rscript 

###############################################################################
###############################################################################

#            PLOTING

###############################################################################
###############################################################################
library(ggplot2)
library(data.table)

v_usage_original <- fread("${original_v_usage}",header =TRUE, data.table = FALSE)

print(v_usage_original[1:4,1:4])

v_usage_tab <- read.csv("${v_usage}",header =TRUE);v_usage_tab[,"X"]=NULL

print(v_usage_tab[1:4,1:4])

# odd columns are of control
chron_v_results <- reshape2::melt(v_usage_original[,seq(1, ncol(v_usage_original), 2) ]) 
chron_v_results[,"variable"] <- gsub("[.]","-", chron_v_results[["variable"]])
chron_v_results[,"method"] <- "control"

# odd columns are of control
chron_v_results_dfn <- reshape2::melt(v_usage_tab[,seq(1, ncol(v_usage_tab), 2) ])
chron_v_results_dfn[,"variable"] <- gsub("[.]","-", chron_v_results_dfn[["variable"]])
chron_v_results_dfn[,"method"] <- "control-DolphinNext"

plot_data <- rbind(chron_v_results, chron_v_results_dfn)

pdf("figure2_panelC_pipeAIRR.pdf",width=7.5,height=7)

ggplot(plot_data, aes(x=as.factor(variable), y=value, fill=method)) + 
  geom_boxplot() +
  #ggpubr::theme_pubclean(base_size = 30) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        legend.position=c(.9,0.95), axis.title = element_text(size = 32)) +
  labs(x = "V gene", y = "Frequency", fill = "") +
  scale_fill_manual(values = c("#61D04F","#DF536B"))
  
dev.off()

"""

}


process crohn_plot_figure2 {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*pdf$/) "pdf/$filename"}
input:
 file basics from g_24_csvFile0_g_32
 file q from g_24_csvFile1_g_32
 file V from g_24_csvFile2_g_32
 file VJ from g_24_csvFile3_g_32

output:
 file "*pdf"  into g_32_outputFilePdf00

script:

"""
#!/usr/bin/env Rscript 

###############################################################################
###############################################################################

#            PLOTING

###############################################################################
###############################################################################

basics=read.csv("basics.csv");basics[,"X"]=NULL
q=read.csv("q.csv");q[,"X"]=NULL
v=read.csv("V.csv");v[,"X"]=NULL
vj=read.csv("VJ.csv");vj[,"X"]=NULL
#colnames(v)=gsub('IGHV','V',gsub('\\.','-',colnames(v)))
#colnames(vj)=gsub('\\.','-',colnames(vj))

pdf("plot.pdf",width=7.5,height=7)
m <- rbind(c(1, 2,3), c(4,4,4),c(5,5,5))
layout(m)
z=basics
par(mar=c(3.5,5,2,1),mgp=c(2.5,0.5,0),cex=0.8)
boxplot(z[z[,"stage"]!='case',1]-2,z[z[,"stage"]=='case',1]-2,
        z[z[,"stage"]!='case',2]-2,z[z[,"stage"]=='case',2]-2,
        z[z[,"stage"]!='case',3]-2,z[z[,"stage"]=='case',3]-2,col=c(3,2,3,2,3,2),xaxt='n',
        ylab='CDR3 AA length',xlab='precentile',las=1,
        at=c(1,2,3.5,4.5,6,7))
axis(side = 1,at = c(1.5,4,6.5),labels=c('10','50','90'),lwd.ticks = T)
legend('topleft',legend = c('control','CD'),fill=c(3,2),box.col = NA)
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
boxplot(z[z[,"stage"]!='case',6],z[z[,"stage"]=='case',6],
        z[z[,"stage"]!='case',5],z[z[,"stage"]=='case',5],
        z[z[,"stage"]!='case',4],z[z[,"stage"]=='case',4],col=c(3,2,3,2,3,2),xaxt='n',
        ylab='V distance from germline',xlab='precentile',las=1,
        at=c(1,2,3.5,4.5,6,7))
axis(side = 1,at = c(1.5,4,6.5),labels=c('10','50','90'),lwd.ticks = T)
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
boxplot(q[q[,"q"]==0&q[,"stage"]=='control',][,"d"],
        q[q[,"q"]==0&q[,"stage"]=='case',][,"d"],
        q[q[,"q"]==1&q[,"stage"]=='control',][,"d"],
        q[q[,"q"]==1&q[,"stage"]=='case',][,"d"],
        q[q[,"q"]==2&q[,"stage"]=='control',][,"d"],
        q[q[,"q"]==2&q[,"stage"]=='case',][,"d"],
        q[q[,"q"]==3&q[,"stage"]=='control',][,"d"],
        q[q[,"q"]==3&q[,"stage"]=='case',][,"d"],
        q[q[,"q"]==4&q[,"stage"]=='control',][,"d"],
        q[q[,"q"]==4&q[,"stage"]=='case',][,"d"],
        col=rep(c(3,2),4),xaxt='n',ylab='diversity',xlab='q'  ,las=1 ,
        at=c(1,2,3.5,4.5,6,7,8.5,9.5,11,12))
axis(side = 1,at = c(1.5,4,6.5,9,11.5),labels=c(0,1,2,3,4),lwd.ticks = T)
mtext('C', side = 3, line = 0.5, adj = 0, cex = 1.1)
par(cex=0.75)
p=v
w=names(p)
p=p[,order(w)]


par(mar=c(6,5,1,1),mgp=c(4.25,0.5,0))
boxplot(p,xaxt='n',col=rep(c(3,2),50),at=sort(c(1+c(0:49)*2.5,2+c(0:49)*2.5)),
        ylab='frequency',xlab='V gene',las=1)
axis(side = 1,at = 2.5*c(0:49)+1.5,labels=names(p)[2*1:50-1],lwd.ticks = T,las=2)
mtext('D', side = 3, line = 1, adj = 0, cex = 1.1)

p=vj

w=names(p)
p=p[,order(w)]
par(mar=c(6,5,1,1))
boxplot(p,xaxt='n',col=rep(c(3,2),50),ylab='frequency',xlab='V & J gene',las=1
        ,at=sort(c(1+c(0:49)*2.5,2+c(0:49)*2.5)))
axis(side = 1,at = 2.5*c(0:49)+1.5,labels=gsub('IGH','',names(p)[2*1:50-1]),lwd.ticks = T,las=2)
mtext('E', side = 3, line = 1, adj = 0, cex=1.1)
dev.off()

"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
