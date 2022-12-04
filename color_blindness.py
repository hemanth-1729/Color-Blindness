# This program takes more than 1 hour to run, might vary depending on the processor
#To run the program under 5mins
#Edit line 230 range to 2936000-2942001 instead of 0 to len(reads)
#because all the reads that match to green or red exons belong to that range

import numpy as np
import pandas as pd
import time
import math
bwt=''
ref=''
reads=[]
# storing Reference Sequence
with open('chrX.fa') as my_file:
    while True:
        ref_line=my_file.readline()
        if not ref_line:
            break
        else:
            ref=ref+ref_line.strip()
        
my_file.close()
ref=ref[5:] # sequence starts after header

#storing Burrows-wheeler transform
with open('chrX_last_col.txt') as my_file:
    while True:
        bwt_line=my_file.readline()
        if not bwt_line:
            break
        else:
            bwt=bwt+bwt_line.strip()
my_file.close()

#storing mapping of indexes in bwt with index in ref
pd_bwt_map=pd.read_csv('chrX_map.txt', sep=" ", header=None)
bwt_map=pd_bwt_map.to_numpy()

#storing the reads
with open('reads') as my_file:
    while True:
        read_line=my_file.readline()
        if not read_line:
            break
        else:
            reads.append(read_line.strip())
my_file.close()


#compute ranks
rank=np.zeros((151100561,4),dtype='int32')
count_A=0
count_C=0
count_G=0
count_T=0
for i in range(len(bwt)):
    if(bwt[i]=='A'):
        count_A+=1
        rank[i]=[count_A,count_C,count_G,count_T]
    elif(bwt[i]=='C'):
        count_C+=1
        rank[i]=[count_A,count_C,count_G,count_T]
    elif(bwt[i]=='G'):
        count_G+=1
        rank[i]=[count_A,count_C,count_G,count_T]
    elif(bwt[i]=='T'):
        count_T+=1
        rank[i]=[count_A,count_C,count_G,count_T]
        
#Get complimentary read        
def get_complimentary(test_read):
    comp=''
    for i in range(len(test_read)):
        if(test_read[i]=='A'):
            comp+='T'
        elif(test_read[i]=='C'):
            comp+='G'
        elif(test_read[i]=='G'):
            comp+='C'
        elif(test_read[i]=='T'):
            comp+='A'
    comp=comp[::-1]
    return comp

#Update bands
def update_band(band,char):
    if(char=='A'):
        low_rank=rank[band[0]][0]
        high_rank=rank[band[1]][0]
    elif(char=='C'):
        low_rank=rank[band[0]][1]
        high_rank=rank[band[1]][1]
    elif(char=='G'):
        low_rank=rank[band[0]][2]
        high_rank=rank[band[1]][2]
    elif(char=='T'):
        low_rank=rank[band[0]][3]
        high_rank=rank[band[1]][3]
    if(low_rank==high_rank and bwt[band[0]]!=char):
        return [0,0]
    else:
        if(char=='A'):
            if(low_rank==0):
                low=low_rank
            else:
                low=low_rank-1
            high=high_rank-1
        elif(char=='C'):
            low=count_A+low_rank-1
            high=count_A+high_rank-1
        elif(char=='G'):
            low=count_A+count_C+low_rank-1
            high=count_A+count_C+high_rank-1
        elif(char=='T'):
            low=count_A+count_C+count_G+low_rank-1
            high=count_A+count_C+count_G+high_rank-1
    return [low,high]

#Get starting band
def get_starting_band(char):
    if(char=='A'):
        band=[0,count_A-1]
    elif(char=='C'):
        band=[count_A,count_A+count_C-1]
    elif(char=='G'):
        band=[count_A+count_C,count_A+count_C+count_G-1]
    elif(char=='T'):
        band=[count_A+count_C+count_G,count_A+count_C+count_G+count_T-1]
    return band
                                
#Get matching info that includes band and offset
def get_matching_info(arr):
    mismatch_count=0
    flag=0
    band=get_starting_band(arr[0])
    #print(band,arr[0])
    for i in range(1,len(arr)):
        if(flag==1):
            band=get_starting_band(arr[i])
            flag=0
        else:
            band=update_band(band,arr[i])
        #print(band,arr[i])

        if(mismatch_count<3):
            if(band[0]==0 and band[1]==0):
                mismatch_count+=1
                flag=1   
            elif(band[1]-band[0]<=3):
                return [band,(i+1-len(arr))] #band,offset
        else:
            return [[0,0],(i+1-len(arr))]#band,offset
    return [band,(i+1-len(arr))]
    
        
#To get starting indexes       
def get_starting_index(band,offset,test_read):
    start_index=[]
    if(band[1]!=0):
        for i in range(band[0],band[1]+1):
            curr_start_index=bwt_map[i][0]+offset
            string_1=ref[curr_start_index:(curr_start_index+len(test_read))]
            mismatch_count=sum(c1!=c2 for c1,c2 in zip(string_1,test_read))
            #print(curr_start_index,len(string_1),len(test_read))
            #print(mismatch_count)
            if(mismatch_count<=2):
                start_index.append(curr_start_index)
        return start_index
                
    return start_index


#To find information of corresponding matched exon        
def find_gene(test_read):
    rev=test_read[::-1]
    band,offset=get_matching_info(rev) 
    start_index=get_starting_index(band,offset,test_read)
    #print(start_index)
    exons=np.zeros((2,6))
    if(len(start_index)==0):
        return exons
    else:
        red_exon=[0,0,0,0,0,0]
        green_exon=[0,0,0,0,0,0]
        for i in range(len(start_index)):
            if((149249757<=start_index[i]<=149249868)or(149249757<=(start_index[i]+len(test_read))<=149249868)): 
                red_exon[0]=1  
            elif((149256127<=start_index[i]<=149256423)or(149256127<=(start_index[i]+len(test_read))<=149256423)):
                red_exon[1]=1
            elif((149258412<=start_index[i]<=149258580)or(149258412<=(start_index[i]+len(test_read))<=149258580)):
                red_exon[2]=1  
            elif((149260048<=start_index[i]<=149260213)or(149260048<=(start_index[i]+len(test_read))<=149260213)):
                red_exon[3]=1   
            elif((149261768<=start_index[i]<=149262007)or(149261768<=(start_index[i]+len(test_read))<=149262007)):
                red_exon[4]=1  
            elif((149264290<=start_index[i]<=149264400)or(149264290<=(start_index[i]+len(test_read))<=149264400)):
                red_exon[5]=1
            elif((149288166<=start_index[i]<=149288277)or(149288166<=(start_index[i]+len(test_read))<=149288277)):
                green_exon[0]=1
            elif((149293258<=start_index[i]<=149293554)or(149293258<=(start_index[i]+len(test_read))<=149293554)):
                green_exon[1]=1
            elif((149295542<=start_index[i]<=149295710)or(149295542<=(start_index[i]+len(test_read))<=149295710)):
                green_exon[2]=1
            elif((149297178<=start_index[i]<=149297343)or(149297178<=(start_index[i]+len(test_read))<=149297343)):
                green_exon[3]=1
            elif((149298898<=start_index[i]<=149299137)or(149298898<=(start_index[i]+len(test_read))<=149299137)):
                green_exon[4]=1
            elif((149301420<=start_index[i]<=149301530)or(149301420<=(start_index[i]+len(test_read))<=149301530)):
                green_exon[5]=1

        for i in range(6):
            if(red_exon[i]==1 and green_exon[i]==1):
                exons[0][i]=0.5
                exons[1][i]=0.5
                return exons
            elif(red_exon[i]==1):
                exons[0][i]=1
                return exons
            elif(green_exon[i]==1):
                exons[1][i]=1
                return exons
        return exons
    
#Intialize exon_count
exon_count=np.zeros((2,6))
begin=time.time()


#Matching reads                
for i in range(len(reads)):
    curr_read=''
    for j in range(len(reads[i])):
        if(reads[i][j]=='N'):
            curr_read+='A'#To replace 'N' with 'A'
        else:
            curr_read+=reads[i][j] 
    exons=find_gene(curr_read)
    if(np.sum(exons)!=0):
        exon_count+=exons
    else:
        comp=get_complimentary(curr_read)
        exons=find_gene(comp)
        exon_count+=exons
    #if(i%10000==0):
        #print(exon_count,i)
end=time.time()

print("---------------------------\n\n")          
print("Time Taken to run program is :",end-begin)   
print("Number of red exons reads are:",exon_count[0][0],exon_count[0][1],
      exon_count[0][2],exon_count[0][3],exon_count[0][4],exon_count[0][5])
print("Number of Green exons reads are:",exon_count[1][0],exon_count[1][1],
      exon_count[1][2],exon_count[1][3],exon_count[1][4],exon_count[1][5]) 
print("\n----------------------------\n") 

#To determine the probability of model based on counts
def generate_prob(model,exon_count):
    prob=1
    for i in range(4):
        a=int(round(exon_count[0][i+1]))
        b=int(round(exon_count[1][i+1]))
        curr_exon_prob=math.comb(a+b,a)*(model[i]**a)*((1-model[i])**b)
        prob=prob*curr_exon_prob
    return prob

#Each list is a model where each element in list indicate probability of a read 
#coming from the corresponding red exon
models=[[1/3,1/3,1/3,1/3],[1/2,1/2,0,0],[1/4,1/4,1/2,1/2],[1/4,1/4,1/4,1/2]] 

prob_list=[]

for i in range(4):
    prob_list.append(generate_prob(models[i],exon_count))

best_prob=max(prob_list)

for i in range(4):
    if(best_prob==prob_list[i]):
        if(i==0):
            print("Person is not colour-blind, Since the Fraction on Red Exons 2,3,4,5 is 50%,50%,50%,50%")
        elif(i==1):
            print("Configuration that lead to colour-blindness is: Fraction on Red Exons 2,3,4,5 is 100%,100%,0%,0%")
        elif(i==2):
            print("Configuration that lead to colour-blindness is: Fraction on Red Exons 2,3,4,5 is 33%,33%,100%,100%")
        else:
            print("Configuration that lead to colour-blindness is: Fraction on Red Exons 2,3,4,5 is 33%,33%,33%,100%")
            
