#!/usr/bin/env python3


def cutSeq(seq, num, cont=0, listSeq=[]):
	if len(seq[cont:])<num:
		return listSeq
	
	cont+=1
	listSeq.append(seq[:num])
	return cutSeq(seq[cont:], num)


def listToFasta(pepList, fastaFormat=None):
	if fastaFormat is None:
		fastaFormat=[]
	[fastaFormat.append(f'>Seq{num}\n{pep}\n') for num, pep in enumerate(pepList, start=1)]
	return fastaFormat


def predictionAMP(fastaFormat, num, data=None, seq=None):
	if data is None:
		data=[]
	if seq is None:
		seq=[]
	from selenium import webdriver
	from selenium.webdriver.chrome.service import Service
	from webdriver_manager.chrome import ChromeDriverManager
	from selenium.webdriver.common.by import By
	from time import sleep
	import pandas as pd
	
	# Define Brave and chromedriver path
	chromedriver = "./chromedriver"
	option = webdriver.ChromeOptions()
	option.binary_location = '/Applications/Brave Browser.app/Contents/MacOS/Brave Browser'
	
	s = Service(chromedriver)
	
	# Create new automated instance of Brave
	browser = webdriver.Chrome(service=s, options=option)
	#browser = webdriver.Chrome(ChromeDriverManager().install(), options=option)
	browser.get("http://www.camp.bicnirrh.res.in/predict/")
	
	# Submitting information to the website
	browser.find_element(By.XPATH, value='//*[@id="frm1"]/p[1]/textarea').send_keys(fastaFormat)
	sleep(2)
	browser.find_element(By.XPATH, value='//*[@id="frm1"]/p[6]/label/input').click()
	sleep(2)
	browser.find_element(By.XPATH, value='//*[@id="frm1"]/p[7]/input[1]').click()
	sleep(5)
	
	# Extracting the data of the table
	table1 = browser.find_element(By.XPATH, value='//table/tbody/tr/td/table[3]/tbody/tr/td/table[1]').text
	table2= browser.find_element(By.XPATH, value='//table/tbody/tr/td/table[3]/tbody/tr/td/table[2]').text
	table3= browser.find_element(By.XPATH, value='//table/tbody/tr/td/table[3]/tbody/tr/td/table[3]').text
	table4= browser.find_element(By.XPATH, value='//table/tbody/tr/td/table[3]/tbody/tr/td/table[4]').text
	
	# Removing unnecessary information
	svm=table1[31:].split('\n')
	rf=table2[31:].split('\n')
	ann=table3[15:].split('\n')
	da=table4[31:].split('\n')
	
	# Transforming data into row
	[data.append(i.split()) for i in svm]
	df1 = pd.DataFrame(data, columns=['ID','SVM', 'SVM_Score'])
	data.clear()
	[data.append(i.split()) for i in rf]
	df2 = pd.DataFrame(data, columns=['ID', 'RF', 'RF_Score'])
	df2.drop('ID', axis=1, inplace=True)
	data.clear()
	[data.append(i.split()) for i in ann]
	df3 = pd.DataFrame(data, columns=['ID', 'ANN'])
	df3.drop('ID', axis=1, inplace=True)
	data.clear()
	[data.append(i.split()) for i in da]
	df4 = pd.DataFrame(data, columns=['ID', 'DA', 'DA_Score'])
	df4.drop('ID', axis=1, inplace=True)
	# Creating the dataframe
	big_df=pd.concat([df1, df2, df3, df4], ignore_index=False, axis=1)
	
	#Seleciona apenas as sequências 
	[seq.append(line[-(num+1):].replace('\n', '')) for line in fastaFormat]
	#Insere as sequências no Dataframe na posição 1
	big_df.insert(1, "Sequence", seq, allow_duplicates=False)
	#Filtra a tabela em duas df_Pep para os que apresentaram predição antimicrobiana e df_NPep para os que não apresentaram
	df_Pep =  big_df.query('SVM=="AMP" & RF=="AMP" & ANN=="AMP" & DA=="AMP"')
	df_Pep.drop('ID', axis=1, inplace=True)
	df_NPep = big_df.query('SVM=="NAMP" or RF=="NAMP" or ANN=="NAMP" or DA=="NAMP"')
	
	import os.path
	#Verifica se o arquivo PepCandidatos existe, caso sim ele adiciona os novos resultados, caso não ele cria o arquivo
	if(os.path.isfile('PepCandidatos.csv')):
		df = pd.read_csv('PepCandidatos.csv')
		df_new = pd.concat([df,df_Pep])
		df_new.to_csv('PepCandidatos.csv', index=False)
	else:
		df_Pep.to_csv('PepCandidatos.csv', index=False)
		
	#Retorna os peptídeos que não apresentaram atividade
	
	return df_NPep['Sequence']

def modPep(listSeq, listPepMod=None):
	if listPepMod is None:
		listPepMod=[]
	import random
	listAmino=['X','X','X','X']
	listSwap=['X', 'X', 'X', 'X']
	
	for seq in listSeq:
		for num, amino in enumerate(seq):
			if amino in listAmino:
				seq=seq.replace(seq[num], random.choice(listSwap),1)
				listPepMod.append(seq)
	return listPepMod

#def(list1, *args):
def modMotif(listSeq, listPepMod=None):
	if listPepMod is None:
		listPepMod=[]
	import random
	import pandas as pd
	listMotif=['XXXXX', 'XXXXX', 'XXXXX', 'XXXXX', 'XXXXX', 'XXXXX', 'XXXXX', 'XXXXX', 'XXXXX', 'XXXXX']
	num=1
	df=pd.DataFrame()
	for seq in listSeq:
		for pos in range(len(seq)-5):
			if seq not in listPepMod:
				listPepMod.append(seq[:pos]+random.choice(listMotif)+seq[pos+5:])
				df=df.append(listPepMod, ignore_index=True)
				
				listPepMod.clear()
				
	[listPepMod.append(f'>Seq{num}\n{pep}\n') for num, pep in enumerate(df[0], start=1)]
	
	return listPepMod

if __name__ == '__main__':
	step1=cutSeq('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX', 20)
	step2=listToFasta(step1)
	step3=predictionAMP(step2, 20)
	print('Approach 1 completed')
	step4=modPep(step3)
	step5=listToFasta(step4)
	step6=predictionAMP(step5, 20)
	print('Approach 2 completed')
	print(step6)
	step7=modMotif(step6)
	print(step7)
	#step8=listToFasta(step7)
	#step9=predictionAMP(step7, 20)
	step8=predictionAMP(step7, 20)
	#print(step7)
	print('Approach 3 completed')
