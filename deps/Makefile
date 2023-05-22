.PHONY: all
all: labw_utils yasim

labw_utils:
	git clone --branch devel https://codeup.aliyun.com/646607f7b60a10415c320130/Bio/labw_utils.git

yasim:
	git clone --branch devel https://codeup.aliyun.com/646607f7b60a10415c320130/Bio/yasim.git

.PHONY: update
update:
	cd labw_utils && git pull && cd ..
	cd yasim && git pull && cd ..
