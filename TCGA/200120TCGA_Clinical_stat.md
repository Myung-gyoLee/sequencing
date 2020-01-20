/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical/02cancer33_txt

## loading data
```r
txtdir = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical"
setwd(txtdir)

cl_files=list.files("02cancer33_txt", pattern="*.txt", full.names=TRUE)


```

## Make column to Dictionary
```python
import os, glob, pandas as pd
txtdir = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical/02cancer33_txt"
txtfiles=glob.glob("%s/*.txt"%txtdir)
txtfiles


heads={}
for tfile in txtfiles:
    tf=open(tfile)
    rline=tf.readline()
    heads[tfile.split("/")[-1].split("_")[-1].split(".")[0]]=rline
    tf.close()

```

##
```python



```

##
```python



```

##
```python



```


##
```python



```

##
```python



```



##
```r

```


##
```r

```


##
```r

```


##
```r

```
