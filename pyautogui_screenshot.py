#-*-coding: utf-8

import pyautogui
import time


## txt 파일 읽어서 gene name 가져오기
in_gene = "MALAT1"

def cap_ccle(in_gene,cluster): 
	import pyautogui
	import time
	## 마우스로 주소 입력 창 클릭
	pyautogui.click(x=449, y=108)
	## double-click address
	pyautogui.click(clicks = 2, x=449, y=108)
	## 주소 입력 
	ad_ccle="page?gene=%s"%(in_gene)
	pyautogui.typewrite(ad_ccle, interval=0.1)
	pyautogui.press("enter")
	time.sleep(5)
	## click down x2 , up x1
	pyautogui.press('down')
	time.sleep(1)
	pyautogui.press('down')
	time.sleep(1)
	pyautogui.press('down')
	time.sleep(1)
	pyautogui.press('up')
	time.sleep(5)
	## screenshot
	im1=pyautogui.screenshot('/home/cytogenbi2/SingleCell/%s_%s.png'%(in_gene,cluster))
	time.sleep(5)
	## click anywhere
	pyautogui.click(x=126, y=637)
	# 안전모드 설정하기, 잘못되었을 경우 탈출구
	pyautogui.PAUSE = 1  
	pyautogui.FAILSAFE = True

if len(in_gene) > 1:
else:break  
	
	


glist=[]
with open("/home/cytogenbi2/SingleCell/CMC009_log2fc2cut.txt") as f:
	for gline in f:
		gene_cluster = "%s:%s"%(gline.split("\t")[0],gline.split("\t")[7])
		glist.append(gene_cluster)
		
glist=glist[5:]
for gline in glist:
	in_gene = gline.split(":")[0]
	cluster = gline.split(":")[1]
	if len(in_gene) > 1:
		cap_ccle(in_gene,cluster)
	else:break
	
from PIL import Image
import glob
import os

inf = glob.glob("/home/cytogenbi2/SingleCell//*.png")
for infile in inf:
	img = Image.open(infile)
	area = (321,449,1728,1061)
	cropped_img = img.crop(area)
	#cropped_img.show()
	cropped_img.save(infile)



def cap_gepia(in_gene,cluster): 
	import pyautogui
	import time
	# 마우스로 gene name 입력 칸 클릭하기
	## 반응형 1868 X 884 
	time.sleep(2)
	#pyautogui.click(x=493, y=563) 
	pyautogui.click(x=541, y=609) 
	#pyautogui.click(x=194, y=525)
	##gene name을 keyboard 입력
	time.sleep(2)  
	pyautogui.typewrite(in_gene, interval=0.25)
	time.sleep(2)
	## click GoPIA!
	## 창 최대화 
	pyautogui.click(x=1410, y=620)
	#pyautogui.click(x=1387, y=565)
	#pyautogui.click(x=747, y=530)
	time.sleep(2)
	## log2 "on"
	pyautogui.click(x=1287, y=777)
	time.sleep(5)
	## page down button
	pyautogui.press('pagedown')
	#pyautogui.dragTo(x=100, y=100)
	## up button 
	pyautogui.press('up')
	time.sleep(5)
	## screenshot
	im1=pyautogui.screenshot('/home/cytogenbi2/SingleCell/gepia_plot/%s_%s.png'%(in_gene,cluster))
	## click address
	pyautogui.click(x=449, y=108)
	pyautogui.click(clicks = 2, x=449, y=108)
	## double-click address
	ad_gepia="http://"
	pyautogui.typewrite(ad_gepia, interval=0.25)
	pyautogui.press("enter")
	## press pop-up "close"
	time.sleep(5)
	pyautogui.click(x=126, y=637)
	# 안전모드 설정하기, 잘못되었을 경우 탈출구
	pyautogui.PAUSE = 1  
	pyautogui.FAILSAFE = True  


glist=[]
with open("/home/cytogenbi2/SingleCell/CMC009_log2fc2cut.txt") as f:
	for gline in f:
		gene_cluster = "%s:%s"%(gline.split("\t")[0],gline.split("\t")[7])
		glist.append(gene_cluster)
		

for gline in glist:
	in_gene = gline.split(":")[0]
	cluster = gline.split(":")[1]
	cap_gepia(in_gene,cluster)


img=Image.open(inf[1])
area = (432,312,1574,969)
cropped_img = img.crop(area)
cropped_img.show()


inf = glob.glob("/home/cytogenbi2/SingleCell/gepia_plot/*.png")
for infile in inf:
	img = Image.open(infile)
	area = (432,314,1572,1375)
	cropped_img = img.crop(area)
	#cropped_img.show()
	cropped_img.save(infile)



'''
#전체 화면 크기
pyautogui.size()

#현재 마우스 위치 확인
pyautogui.position()


# 모니터 해상도 가져오기
width, height = pyautogui.size()  
print('width={0}, height={1}'.format(width, height))  


# 마우스 위치 가져오기
x, y = pyautogui.position()  
print('x={0}, y={1}'.format(x, y))  

# 안전모드 설정하기, 잘못되었을 경우 탈출구
pyautogui.PAUSE = 1  
pyautogui.FAILSAFE = True  
'''

