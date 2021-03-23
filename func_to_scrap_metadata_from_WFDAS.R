#func to scrap metadata from WFDAS
require(rvest)
scrapWFDASmd<-function(site,gacc,state,group){
	#get the url
	#require(rvest)
	url<-paste0('https://www.wfas.net/nfmd/include/site_page.php?site=',
				site,
				'&gacc=',gacc,
				'&state=',state,
				'&grup=',group)
				
	url<-gsub(" ", "%20", url)
	#scrap for the md
	page<-read_html(url)
	site_parts<-html_nodes(page,'table')
	table<-html_children(site_parts)
	coords<-html_text(table[[4]],trim = T)
	elev<-html_text(table[[5]])
	elev<-do.call(rbind,strsplit(elev, "\\n"))[,2]	
	elev<-gsub("[[:space:]]", "",elev)
	elev<-as.numeric(gsub(",", "",elev))
	coords<-do.call(rbind,strsplit(coords, "\\n"))
	coords<-do.call(rbind,strsplit(coords[1,2], "x",fixed = T))
	if(length(coords[1,])<2){
		lat<-NA
		lon<-NA
	}else{
		
		lat<-do.call(rbind,strsplit(coords[1,1], "-",fixed = T))	
		lon<-do.call(rbind,strsplit(coords[1,2], "-",fixed = T))
		lat<-as.numeric(gsub("[[:space:]]", "",lat[1,1]))+ as.numeric(gsub("[[:space:]]", "",lat[1,2]))/60 + as.numeric(gsub("[[:space:]]", "",lat[1,3]))/60
		
		lon<-as.numeric(gsub("[[:space:]]", "",lon[1,1]))+ as.numeric(gsub("[[:space:]]", "",lon[1,2]))/60 + as.numeric(gsub("[[:space:]]", "",lon[1,3]))/60
	}

	data.frame(lat=lat,lon=-1*lon,elev=elev)
	
}
















