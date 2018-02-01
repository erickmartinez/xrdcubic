var nPeaks = 0;

function precisionRound(number, precision) {
    var factor = Math.pow(10, precision);
    return Math.round(number * factor) / factor;
}

var getMillerIndices = function(n) {
    var hashIndices = new Map();
    var sorted = [];
    var count = 0;
    var d0 = 1;
    for (var i=1;i<=n;++i) {
        for (var j=0;j<n;++j) {
            for (var k=0;k<n;++k){
                var keyarr = [i,j,k];
                keyarr = keyarr.sort(function(a, b){return b-a})
                var key = keyarr.join('');
                if (!(key in hashIndices)) {
                    var mIndex = new Object();
                    mIndex.fcc = isFCC(i,j,k);
                    mIndex.bcc = isBCC(i,j,k);
                    mIndex.d = 1.0/Math.sqrt(i*i+j*j+k*k);
                    if (count == 0) {
                        d0 = mIndex.d;
                    }
                    mIndex.ratio = d0/mIndex.d;
                    mIndex.index = key;
                    hashIndices[key] = mIndex;
                    ++count;
                }
            }
        }
    }
    for (midx in hashIndices) {
        sorted.push(hashIndices[midx]);
    }

    sorted.sort(function(a,b){
        return b.d - a.d;
    });
    return sorted;
};

var getFCCIndices = function(n) {
    var hashIndices = new Map();
    var sorted = [];
    var count = 0;
    var d0 = 0;
    for (var i=1;i<=n;++i) {
        for (var j=0;j<n;++j) {
            for (var k=0;k<n;++k){
                var keyarr = [i,j,k];
                keyarr = keyarr.sort(function(a, b){return b-a})
                var key = keyarr.join('');
                if (!(key in hashIndices)) {
                    var mIndex = new Object();
                    if (isFCC(i,j,k)) {
                        mIndex.d = 1.0/Math.sqrt(i*i+j*j+k*k);
                        if (count == 0) {
                            d0 = mIndex.d;
                        }
                        mIndex.ratio = d0/mIndex.d;
                        mIndex.index = key;
                        hashIndices[key] = mIndex;
                        ++count;
                    }                    
                }
            }
        }
    }
    for (midx in hashIndices) {
        sorted.push(hashIndices[midx]);
    }

    sorted.sort(function(a,b){
        return b.d - a.d;
    });
    return sorted;
}

var getBCCIndices = function(n) {
    var hashIndices = new Map();
    var sorted = [];
    var count = 0;
    var d0 = 0;
    for (var i=1;i<=n;++i) {
        for (var j=0;j<n;++j) {
            for (var k=0;k<n;++k){
                var keyarr = [i,j,k];
                keyarr = keyarr.sort(function(a, b){return b-a})
                var key = keyarr.join('');
                if (!(key in hashIndices)) {
                    var mIndex = new Object();
                    if (isBCC(i,j,k)) {
                        mIndex.d = 1.0/Math.sqrt(i*i+j*j+k*k);
                        if (count == 0) {
                            d0 = mIndex.d;
                        }
                        mIndex.ratio = d0/mIndex.d;
                        hashIndices[key] = mIndex;
                        ++count;
                    }                    
                }
            }
        }
    }
    for (midx in hashIndices) {
        sorted.push(hashIndices[midx]);
    }

    sorted.sort(function(a,b){
        return b.d - a.d;
    });
    return sorted;
}


var checkFCC = function(peaks,lam) {
    var n = peaks.length;
    var planeDistances = getPlaneDistances(peaks,lam);
    var distanceRatios = getDistanceRatios(planeDistances);
    var milleridx = getFCCIndices(n);
    //console.log(['fcc',milleridx,distanceRatios]);
    var results = [];
    var count = 0;
    var d0 = 1;
    for (var mkey in milleridx) {
        midx = milleridx[mkey];
        dRatio = distanceRatios[count];
        mRatio = midx.ratio;
        if (precisionRound(dRatio,2) == precisionRound(mRatio,2)) {
            oIdx = new Object();
            oIdx.index = midx.index;
            oIdx.d = planeDistances[count];
            oIdx.ratio = dRatio;
            oIdx.angle = peaks[count];
            results.push(oIdx);
            ++count;
        }  
    }
    results.sort(function(a,b){return a.ratio - b.ratio});
    if (results.length > 1) {
        return results;
    }
    return -1;
};

var checkSimple = function(peaks,lam) {
    var n = peaks.length;
    var planeDistances = getPlaneDistances(peaks,lam);
    var distanceRatios = getDistanceRatios(planeDistances);
    var milleridx = getMillerIndices(n);
    //console.log(['simple',milleridx,distanceRatios]);
    var results = [];
    var count = 0;
    var d0 = 1;
    for (var mkey in milleridx) {
        midx = milleridx[mkey];
        dRatio = distanceRatios[count];
        mRatio = midx.ratio;
        if (precisionRound(dRatio,2) == precisionRound(mRatio,2)) {
            oIdx = new Object();
            oIdx.index = midx.index;
            oIdx.d = planeDistances[count];
            oIdx.ratio = dRatio;
            oIdx.angle = peaks[count];
            results.push(oIdx);
            ++count;
        }  
    }
    results.sort(function(a,b){return a.ratio - b.ratio});
    if (results.length > 1) {
        return results;
    }
    return -1;
};

var checkBCC = function(peaks,lam) {
    var n = peaks.length;
    var planeDistances = getPlaneDistances(peaks,lam);
    var distanceRatios = getDistanceRatios(planeDistances);
    var milleridx = getBCCIndices(n);
    //console.log(['bcc',milleridx,distanceRatios]);
    var results = [];
    var count = 0;
    var d0 = 1;
    for (var mkey in milleridx) {
        midx = milleridx[mkey];
        dRatio = distanceRatios[count];
        mRatio = midx.ratio;
        if (precisionRound(dRatio,2) == precisionRound(mRatio,2)) {
            oIdx = new Object();
            oIdx.index = midx.index;
            oIdx.d = planeDistances[count];
            oIdx.ratio = dRatio;
            oIdx.angle = peaks[count];
            results.push(oIdx);
            ++count;
        }  
    }
    results.sort(function(a,b){return a.ratio - b.ratio});
    if (results.length > 1) {
        return results;
    }
    return -1;
};

var isFCC = function(h,k,l){
    if (h%2 == k%2 && k%2 == l%2){
        return true;
    }
    return false;
}

var isBCC = function(h,k,l){
    var sum = h + k + l;
    if (sum%2==0){
        return true;
    }
    return false;
}


var addPeak = function(){
    var $lastARow = $("#tPeaks tr:last");
    var lastId = $lastARow.attr('id');
    var numberId =  (nPeaks == 0) ? 1 : parseInt(lastId.match(/\d+$/)) + 1;
    var rowId = 'row-'+numberId;
    var inputId = 'peak-' + numberId;
    var $row = $("<tr></tr>",{id:rowId});
    $lastARow = $row.insertAfter($lastARow);
    var $td1 = $('<td>Peak #'+numberId+'</td>');
    $lastARow.append($td1);
    var $td2 = $('<td></td>');
    var $inputPos = $('<input/>',{id:inputId});
    $inputPos.attr({'type':"number",
                        'class':"peak-input",
						"step":0.1, 
                        "min":0, "max":120});
    $td2.append($inputPos);
    $lastARow.append($td2);
    var $td3 = $('<td></td>');
    var $delBtn = $('<input/>');
    $delBtn.attr({'value':'delete',type:'button',onclick:'deletePeak("'+rowId+'")'});
    $td3.append($delBtn);
    $lastARow.append($td3);
    nPeaks++;
}

var deletePeak = function(id) {
    if (nPeaks > 0) {
        $('#'+id).remove();
        nPeaks--;
    }
};

$("#addBtn").click(function(){
    addPeak();
});

var getPlaneDistances = function(angles,lam) {
    var distances = [];
    for (key in angles) {
        theta = 0.5*angles[key];
        rad = theta*Math.PI/180.0;
        d = lam/(2.0*Math.sin(rad));
        distances.push(d);
    }
    
    distances.sort(function(a,b){return b-a});
    return distances;
};

var getDistanceRatios = function(distances) {
    d0 = Math.max.apply(null,distances);
    ratios = [];
    for (var i=0; i<distances.length;++i) {
        ratios.push(d0/distances[i]);
    }
    return ratios;
};

$("#indexPeaksBtn").click(function(){
    var peaks = [];
    var lam = $('#lambda').val();
    $(".peak-input").each(function(key,val){
        peaks.push(val.value);
    });
    var peakIndices = checkSimple(peaks,lam);
    if (peakIndices == -1) {
        peakIndices = checkFCC(peaks,lam);
        if (peakIndices == -1) {
            peakIndices = checkBCC(peaks,lam);
            if (peakIndices == -1) {
                printResults([],'Not a cubic system!');
            } else {
                printResults(peakIndices,'BCC');
            }
        } else {
            printResults(peakIndices,'FCC');
        }
    } else {
        printResults(peakIndices,'simple');
    }

});

var getLatticeConstant = function(h,k,l,d) {
    return Math.abs(d)*Math.sqrt(h*h+k*k+l*l);
};

var printResults = function(results,lattice='simple') {
    $("#result-lattice-system").html(lattice);
    $("#xrd-results-tbl tr:gt(0)").remove();
    var p0 = results[0];
    var midx = p0.index;
    var h = parseInt(midx.substring(0,1));
    var k = parseInt(midx.substring(1,2));
    var l = parseInt(midx.substring(2,3));
    var d = p0.d;
    var a = precisionRound(getLatticeConstant(h,k,l,d),3);
    $("#result-lattice-constant").html(a+" &angst;");
    console.log(results);
    for (p in results) {
        var peak = results[p];
        var $lastARow = $("#xrd-results-tbl tr:last");
        var $row = $("<tr></tr>");
        $tr = $row.insertAfter($lastARow);
        var $td1 = $('<td>'+peak.angle+'</td>');
        $tr.append($td1);
        var $td2 = $('<td>'+peak.index+'</td>');
        $tr.append($td2);
        var distance = precisionRound(peak.d,3);
        var $td3 = $('<td>'+distance+'</td>');
        $tr.append($td3);
        var ratio = precisionRound(peak.ratio,3);
        var $td4 = $('<td>'+ratio+'</td>');
        $tr.append($td4);
    }    
}