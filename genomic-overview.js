/*
 * Copyright (c) 2015 Memorial Sloan-Kettering Cancer Center.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF MERCHANTABILITY OR FITNESS
 * FOR A PARTICULAR PURPOSE. The software and documentation provided hereunder
 * is on an "as is" basis, and Memorial Sloan-Kettering Cancer Center has no
 * obligations to provide maintenance, support, updates, enhancements or
 * modifications. In no event shall Memorial Sloan-Kettering Cancer Center be
 * liable to any party for direct, indirect, special, incidental or
 * consequential damages, including lost profits, arising out of the use of this
 * software and its documentation, even if Memorial Sloan-Kettering Cancer
 * Center has been advised of the possibility of such damage.
 */

/*
 * This file is part of cBioPortal.
 *
 * cBioPortal is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

function GenomicOverviewConfig(nRows,width) {
    this.nRows = nRows;
    this.canvasWidth = width;
    this.wideLeftText = 25;
    this.wideRightText = 35;
    this.GenomeWidth = this.canvasWidth-this.wideLeftText-this.wideRightText;
    this.pixelsPerBinMut = 3;
    this.rowHeight = 20;
    this.rowMargin = 5;
    this.ticHeight = -0.3;
    this.cnTh = [0.2,1.5];
    this.cnLengthTh = 50000;
}
GenomicOverviewConfig.prototype = {
    getCnColor: function(cnValue) {
        if (cnValue>=this.cnTh[1])
            return "#f00";
        if (cnValue<=-this.cnTh[1])
            return "#00f";
        var c = Math.round(255*(this.cnTh[1]-Math.abs(cnValue))/(this.cnTh[1]-this.cnTh[0]));
        if (cnValue<0)
            return "rgb("+c+","+c+",255)";
        else
            return "rgb(255,"+c+","+c+")";
    },
    canvasHeight: function() {
        return 2*this.rowMargin+this.ticHeight+this.nRows*(this.rowHeight+this.rowMargin);
    },
    yRow: function(row) {
        return 2*this.rowMargin+this.ticHeight+row*(this.rowHeight+this.rowMargin);
    },
    xRightText: function() {
        return this.wideLeftText + this.GenomeWidth+5;
    }
};

function getChmEndsPerc(chms, total) {
    var ends = [0];
    for (var i=1; i<chms.length; i++) {
        ends.push(ends[i-1]+chms[i]/total);
    }
    return ends;
}

/**
 * storing chromesome length info
 */
function ChmInfo() {
    this.hg19 = [0,249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566];
    this.total = 3095677412;
    this.perc = getChmEndsPerc(this.hg19,this.total);
}
ChmInfo.prototype = {
    loc2perc : function(chm,loc) {
        return this.perc[chm-1] + loc/this.total;
    },
    loc2xpixil : function(chm,loc,goConfig) {
        return this.loc2perc(chm,loc)*goConfig.GenomeWidth+goConfig.wideLeftText;
    },
    perc2loc : function(xPerc,startChm) {
        var chm;
        if (!startChm) {//binary search
            var low = 1, high = this.hg19.length-1, i;
            while (low <= high) {
                i = Math.floor((low + high) / 2); 
                if (this.perc[i] >= xPerc)  {high = i - 1;}
                else  {low = i + 1;}
            }
            chm = low;
        } else {//linear search
            var i;
            for (i=startChm; i<this.hg19.length; i++) {
                if (xPerc<=this.perc[i]) break;
            }
            chm = i;
        }
        var loc = Math.round(this.total*(xPerc-this.perc[chm-1]));
        return [chm,loc];
    },
    xpixil2loc : function(goConfig,x,startChm) {
        var xPerc = (x-goConfig.wideLeftText)/goConfig.GenomeWidth;
        return this.perc2loc(xPerc,startChm);
    },
    middle : function(chm, goConfig) {
        var loc = this.hg19[chm]/2;
        return this.loc2xpixil(chm,loc,goConfig);
    },
    chmName : function(chm) {
        if (chm === 23) return "X";
        if (chm === 24) return "Y";
        return chm;
    }
};

function plotChromosomes(config,chmInfo,elementId) {
    var yRuler = config.rowMargin+config.ticHeight;
    drawLine1(config.wideLeftText,yRuler,config.wideLeftText+config.GenomeWidth,yRuler,'#000',1,elementId);
    // ticks & texts
    for (var i=1; i<chmInfo.hg19.length; i++) {
        var xt = chmInfo.loc2xpixil(i,0,config);
        drawLine2(xt,yRuler,xt,config.rowMargin,'#000',1,elementId);
        
        var m = chmInfo.middle(i,config);
        text(m,yRuler+0.1,chmInfo.chmName(i),elementId);
    }
    drawLine2(config.wideLeftText+config.GenomeWidth,yRuler,config.wideLeftText+config.GenomeWidth,config.rowMargin,'#000',1,elementId);
}

function drawLine1(x11, y11, x22, y22, cl, width1,elementId) {

	var trace1 = {};
	var layout = {
		xaxis: {
		    showticklabels: false,
		    autotick: false,
		    showgrid: false,
		    zeroline: false,
		},
		yaxis: {
		    range: [0, 5],
		    showticklabels: false,
		    autotick: false,
		    showgrid: false,
		    zeroline: false

		},
		width:1200,
		height:500,
		showlegend: false,
		hovermode:false,
		shapes: [
		{
			type: 'line',
			x0: x11,
			y0: y11,
			x1: x22,
			y1: y22,
			line: {
				color: cl,
				width: width1,
			}
		},
		]
	};
	var data = [trace1];
	Plotly.plot(elementId, data, layout);
	return;  
}

function drawLine2(x11, y11, x22, y22, cl, width1,elementId) {
	myDiv.data.push({
		x:[x11,x22],
		y:[y11,y22],
		mode:'lines',
		line: {
			color: cl,
			width: width1,
	  	}
	});  
	Plotly.redraw(elementId);
	return;  
}

function text(x11, y11,number,elementId) {
	myDiv.data.push({
		x:[x11],
		y:[y11],
	  	mode: 'lines+text',
		text: [number],
	  	textposition: 'top',
	  	textfont: {
		  	family: 'sans serif',
			size: 11,
			color: '#ff7f0e'
	  	}

	});  
	Plotly.redraw(elementId);
	return;  
}
