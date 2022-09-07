// Use custom colours for highchart
Highcharts.setOptions({
  colors: [
      '#4572A7', 
      '#AA4643',
      '#444444', 
      '#80699B', 
      '#3D96AE', 
      '#DB843D', 
      '#92A8CD', 
      '#A47D7C', 
      '#9fb876',
      '#a68038']
    });


/* BOXPLOT */

/* CALCULATE BOXPLOT VALUES */
function getBoxValues(values, counts) {
  var boxData = {},
    min = Math.min.apply(Math, values),
    max = Math.max.apply(Math, values),
    q1 = getPercentile(values, 25),
    median = getPercentile(values, 50),
    q3 = getPercentile(values, 75),
    iqr = q3 - q1,
    lowerFence = q1 - (iqr * 1.5),
    upperFence = q3 + (iqr * 1.5);
	
  for (var i = 0; i < values.length; i++) {
    if (values[i] < lowerFence || values[i] > upperFence) {
      outliers.push([values[i], counts[i]]);
    }
  }

  boxData.values = {};
  boxData.values.counts = counts;
  boxData.values.low = lowerFence;
  boxData.values.q1 = q1;
  boxData.values.median = median;
  boxData.values.q3 = q3;
  boxData.values.high = upperFence;
  boxData.outliers = outliers;

  return boxData;
}

// Percentiles from an array
function getPercentile(data, percentile) {
  data.sort(numSort);
  var index = (percentile / 100) * data.length;
  var result;
  if (Math.floor(index) == index) {
    result = (data[(index - 1)] + data[index]) / 2;
  } else {
    result = data[Math.floor(index)];
  }
  return result;
}

// Function to find mean
function findMean(numbers) {
    var total = 0, i;
    for (i = 0; i < numbers.length; i += 1) {
        total += numbers[i];
    }
    return total / numbers.length;
}

// Mean of an array of numbers
function mean(data) {
  var len = data.length;
  var sum = 0;
  for (var i = 0; i < len; i++) {
    sum += parseFloat(data[i]);
  }
  return (sum / len);
}

// Sort
function numSort(a, b) {
  return a - b;
}

// Generate categories
function alphaCats(n) {
  var alph = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'];
  var cats = [];
  for (var i = 0; i < n; i++) {
    if (i < alph.length) {
      cats.push(alph[i]);
    } else {
      var rep = Math.ceil((i / (alph.length - 1)));
      var c = (Math.ceil(i / (rep - 1)) - alph.length);
      var cat = [];
      for (var j = 0; j < rep; j++) {
        cat.push(alph[c]);
      }
      cats.push(cat);
    }
  }
  return cats;
}

/* CALCULATE BOXPLOT VALUES */

function createBoxplot(dataBackground, dataSample, deletion_id) {
	var backgroundValues = dataBackground.map(function(d) {
	return d.ScaledProportion
	});

	var backgroundCounts = dataBackground.map(function(d) {
	return d.Count
	});


	// Testing and temporary way to get ID values (for X axis), and also max (for Y axis)

	let all_values = dataBackground.map(function(d) { return d.ScaledProportion; })
	let max_value = Math.max.apply(Math, all_values)
	let mean_value = findMean(all_values)

	numberOfBoxes = 2,
	numberOfPoints = 100,
	cats = alphaCats(numberOfBoxes), // if needed, atm using background and my sample as cat names
	outliers = [];
	var boxData = [], meanData = [];

	var boxValues = getBoxValues(backgroundValues, backgroundCounts);
	boxData.push(boxValues.ScaledProportion);
	meanData.push([0, mean(backgroundValues)]);
    sample_value = dataSample.ScaledProportion


//Work out max proportion from sample and background for limit 
      const max_s_back = Math.max(...all_values)
      const max_s_plot = Math.max(max_s_back,sample_value)
      const max_f_plot = Math.max(max_s_plot, 1.5)

		Highcharts.chart('plot_box', {
		chart: {
			type: 'boxplot'
		},
		plotOptions: {
			column: {
				colorByPoint: true
			}
		},
    title: {
        text: gene_id + ' ' + deletion_id,
        style: { "fontSize": "14px" }
    },
		legend: {
			enabled: false
		},
		xAxis: {
			title: {
				text: null
			},
		},
		yAxis: {
			title: {
				text: 'ScaledPropotion'
			},
      type: 'logarithmic',
      min: 0.0001,
      max: max_f_plot,
			plotLines: [
        {
          value: 0.001,
          color: '#db24a1', // magenta
          width: 2,
          label: {
            text: 'Low confidence 0.001',
            align: 'center',
            style: {
              color: '#db24a1'
            }
          }
        },
          {
          value: 0.01,
          color: '#ffa500', // orange
          width: 2,
          label: {
            text: 'Medium confidence 0.01',
            align: 'center',
            style: {
              color: '#ffa500'
            }
          }
        },
        {
          value: 0.1,
          color: '#800080', // purple
          width: 4,
          label: {
            text: 'High confidence 0.1',
            align: 'center',
            style: {
                color: '#800080'
            }
          }
        }
      ]
		},
		series: [
			
			{
				name: 'Analysed sample',
				color: Highcharts.getOptions().colors[0],
				type: 'scatter',
				data:  [sample_value],
				zindex: 5,
				marker: {
            radius: 8,
						fillColor: '#10AAEF',
						lineWidth: 3,
						lineColor: '#112591'
				},
				tooltip: {
						pointFormat: 'ScaledProportion: {point.y}'
				}
			}
		],
    credits: {
      enabled: false
    },
    exporting: {
      enabled: false
    }
  });
}
