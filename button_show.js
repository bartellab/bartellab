// JavaScript Document

// handles display of rollover images on = bartel home =

/*
This script and many more are available free online at
The JavaScript Source!! http://javascript.internet.com */

//also need images arrays for the buttons:  one array of "off", one of "on" (5 each)
	buttons_off = new Array;
	buttons_on = new Array;
	
	//array containing preView window images, 480x160
	preViews = new Array;

	//bartel lab home
	buttons_off[0] = new Image;  buttons_off[0].src = "img/miRNA_manu.jpg";  
	buttons_on[0] = new Image;  buttons_on[0].src = "img/processing_summary.jpg";
	
	//Research summary
	preViews[0] = new Image; preViews[0].src = "img/target_summary.jpg";
	
	
	//Lab members
	buttons_off[1] = new Image;  buttons_off[1].src = "img/miRNA_manu.jpg";
	buttons_on[1] = new Image;  buttons_on[1].src = "img/processing_members.jpg";
	preViews[1] = new Image; preViews[1].src = "img/target_members.jpg";
	
    //Alumni
	buttons_off[2] = new Image;  buttons_off[2].src = "img/miRNA_manu.jpg";
	buttons_on[2] = new Image;  buttons_on[2].src = "img/processing_alumni.jpg";
	preViews[2] = new Image; preViews[2].src = "img/target_alumni.jpg";
	
	
	//Publication
	buttons_off[3] = new Image;  buttons_off[3].src = "img/miRNA_manu.jpg";
	buttons_on[3] = new Image;  buttons_on[3].src = "img/processing_publication.jpg";
	preViews[3] = new Image; preViews[3].src = "img/target_publications2.jpg";
	
	
	//Protocols
	buttons_off[4] = new Image;  buttons_off[4].src = "img/miRNA_manu.jpg";
	buttons_on[4] = new Image;  buttons_on[4].src = "img/processing_protocol.jpg";
	preViews[4] = new Image; preViews[4].src = "img/target_protocols.jpg"; 
	
	//Software
	buttons_off[5] = new Image;  buttons_off[5].src = "img/miRNA_manu.jpg";
	buttons_on[5] = new Image;  buttons_on[5].src = "img/processing_softwares.jpg";
	preViews[5] = new Image; preViews[5].src = "img/target_software.jpg";
	
	//RNA&others
	buttons_off[6] = new Image;  buttons_off[6].src = "img/miRNA_manu.jpg";
	buttons_on[6] = new Image;  buttons_on[6].src = "img/processing_targetscan.jpg";
	preViews[6] = new Image; preViews[6].src = "img/target_targetscan.jpg"; 
	
	preViews[7] = new Image; preViews[7].src = "img/target_nothing.jpg";
	function on(imgName, imgNum)  { 
	    //document[imgName].alt="changed"
		//document[imgName].src = "img/croped_summary.jpg";
		document[imgName].src = buttons_on[imgNum].src;
	}
	
	function off(imgName, imgNum) { 
		document[imgName].src = buttons_off[imgNum].src; 
	}
 
	function up(imgName, imgNum) { 
		document[imgName].src = preViews[imgNum].src; 
	}
	
	//display info about lab member onMouseOver in div called 'div_description'
	//string memInfo is new HTML to write
	//props to DUSTIN DIAZ
	function showInfo(memInfo) {
  		/*	adopted from following addElement code
			var ni = document.getElementById('myDiv');
			var numi = document.getElementById('theValue');
		  	var num = (document.getElementById('theValue').value -1)+ 2;
			numi.value = num;
		  	var newdiv = document.createElement('div');
			var divIdName = 'my'+num+'Div';
			newdiv.setAttribute('id',divIdName);
			newdiv.innerHTML = 'Element Number '+num+' has been added! <a href='#' onclick='removeElement('+divIdName+')'>Remove the div "'+divIdName+'"</a>';
			ni.appendChild(newdiv);*/
		var descriptionS = document.getElementById('div_description');
		descriptionS.innerHTML = memInfo;
	}
	

	/* to launch a new window of specified dimensions */
	//source:  www.spoono.com/flash
	function fixedSize(url,name,features) {
		//This launches a new window and then
		//focuses it if window.focus() is supported.
		win = window.open(url,name,features);
	}