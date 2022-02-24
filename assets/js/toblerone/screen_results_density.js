/* DENSITY PLOTS */

// Gaussian Kernel density estimation - https://www.highcharts.com/blog/tutorials/data-science-and-highcharts-kernel-density-estimation/
// function GaussKDE(xi, x) {
//   return (1 / Math.sqrt(2 * Math.PI)) * Math.exp(Math.pow(xi - x, 2) / -2);
// }

// New Gaussian KDE function with bandwidth
function GaussKDE(xi, x, bw = 0.35) {
  return 1/Math.sqrt(2*Math.PI)/bw*Math.exp(-(x-xi)*(x-xi)/2/bw/bw)
}

const precision = 0.0000001 // Below that threshold, don't display KDE dots

function getSamplesThresholdData(data_thresholds, sample_value) {
  data_thresholds_arr = data_thresholds.map(function(d) { return parseFloat(d.Threshold) })

  arr_vals = Object.keys(data_thresholds_arr).filter(el => sample_value >= data_thresholds_arr[el]);
  arr_vals = arr_vals.map(Number)
  max_val = Math.max(...arr_vals)

  threshold_data = data_thresholds[max_val]
  
  return threshold_data;
}

function getMinMax(value) {
return [Math.min(...value), Math.max(...value)]
}

function createRangeArr(start, end, step = 1) {
const len = Math.floor((end - start) / step) + 1
return Array(len).fill().map((_, idx) => start + (idx * step))
}

function createDensityPlot(data_background, sample_value, data_thresholds, data_type, deletion_id) {
  gaussian_points = 10 // Number of points a gaussian curve will be divided to

  // Display values either as scaled proportion or z-scored based on the selected data type
  if (data_type == "proportion") {
    let arr_all_values = [] // Here goes all values (background + sample) which will be used to calculate a maximum value for finding endpoint for the plot

    sample_value_for_plot = sample_value.ScaledProportion

    data_background = data_background.map(function(d) {
      arr_all_values.push(d.ScaledProportion)
      return d.ScaledProportion
    })

    let arr_all_values_min = 0
    let arr_all_values_max = Math.max(...arr_all_values)
    if (arr_all_values_max == 0) { arr_all_values_max = 0.01 } // If max is also 0, then make max 0.01 to be able to plot anything

    if (arr_all_values_max < sample_value_for_plot) { arr_all_values_max = sample_value_for_plot }
    if (arr_all_values_min > sample_value_for_plot) { arr_all_values_min = sample_value_for_plot }
    
    let gbw = (arr_all_values_max > 0.1 ? 0.02 : 0.0005) // Two gaussian bw will be divided for max val > 0.1 and below it (smaller)
    
    var range = createRangeArr(arr_all_values_min, arr_all_values_max, arr_all_values_max/gaussian_points) // For ScaledProportions plot, calculate kernels for each 0.1
    var gaussian_bw = gbw
  }
  else if (data_type == "zscore") {
    let arr_all_values = [] // Here goes all values (background + sample) which will be used to calculate a maximum value for finding endpoint for the plot

    sample_value_for_plot = sample_value.Zscore

    data_background = data_background.map(function(d) { 
      arr_all_values.push(d.Zscore)

      return d.Zscore
    })

    let arr_all_values_minmax = getMinMax(arr_all_values)  // Take the minimum value of all values in array
    let arr_all_values_min = arr_all_values_minmax[0], arr_all_values_max = arr_all_values_minmax[1]

    // if sample is outside cohort bounds, need to make that the margins of plot
    if (arr_all_values_max < sample_value_for_plot) { arr_all_values_max = sample_value_for_plot }
    if (arr_all_values_min > sample_value_for_plot) { arr_all_values_min = sample_value_for_plot }
    if (arr_all_values_max == 0) { arr_all_values_max = 1 } // If max is also 0, then make max 1 to show a plot

    let startPoint = arr_all_values_min, endPoint = Math.ceil(arr_all_values_max)

    var range = createRangeArr(startPoint, endPoint, endPoint/gaussian_points)

    var gaussian_bw = 0.2
  }

  let sample_threshold = getSamplesThresholdData(data_thresholds, sample_value.ScaledProportion) // Match sample value with threshold data to find the name and colour of the threshold it exceeds - use only ScaledProportion

  let dataSource = data_background;
  let xiData = [];
  let kernel = [];
  let animationDuration = 1500;

  range.forEach((i) => {
    xiData[i] = i;
  })

  let data = []; // data array for plot
  let data_N = dataSource.length;

  for (i = 0; i < range.length; i++) {
    let temp = 0;
    for (j = 0; j < data_N; j++) {
      temp = temp + GaussKDE(range[i], dataSource[j], gaussian_bw);
    }

    let kde_val = (1 / data_N) * temp
    if (kde_val > precision) {
      data.push([range[i], kde_val]);
    }
    else {
      data.push([range[i], null]);
    }
  }

  for (i = 0; i < range.length; i++) {
    let temp = 0;
    kernel.push([]);

    kernel[i].push(new Array(data_N));

    for (j = 0; j < data_N; j++) {
        temp = temp + GaussKDE(range[i], dataSource[j], gaussian_bw);
        kernel[i][j] = GaussKDE(range[i], dataSource[j], gaussian_bw);
    }
  }

  // Add BACKGROUND data to series
  let data_scatter = [];

  data_background.forEach((value, i) => {
      data_scatter.push({
              x: value, // X axis value is the zscore
              y: 0, // Display bottom on the Y axis
      })
  })

  // Create the density plot itself
  highcharts_density_plot = Highcharts.chart("plot_density", {
    chart: {
      type: "spline",
      animation: true
    },
    title: {
        text: "Results for " + gene_id + ' ' + deletion_id,
        style: { "fontSize": "14px" }

    },
    xAxis: {
        plotLines: [{
            value: sample_value_for_plot,
            width: 3,
            color: sample_threshold.Colour,
            dashStyle: 'dash',
            zIndex: 20,
            label: {
                text: (sample_threshold.Threshold == -1 ? 'Sample' : sample_threshold.Label), // Don't show label when there is no reference (threshold is -1)
                align: 'left',
                style: {
                    fontSize: '18px',
                    color: sample_threshold.Colour,
                    zIndex: 20
                }
            }
        }]
    },
    yAxis: {
      title: { text: null }
    },
    tooltip: {
      valueDecimals: 6
    },
    plotOptions: {
        series: {
          marker: {
              enabled: true,
              radius: 4,
              fillColor: "#111111"
          },
          dashStyle: "shortdot",
          color: "#5aa1cb",
          pointStart: range[0],
          animation: {
              duration: animationDuration
          }
        }
    },
    series: [
        {
          name: "Cohort Density",
          data: data,
          connectNulls: true,
          dashStyle: "solid",
          lineWidth: 3,
          color: "#5aa1cb"
        },
        {
          type: 'scatter',
          name: 'Cohort data points',
          marker: {
            enabled: true,
            radius: 10,
            symbol: "triangle",
            fillColor: "#5aa1cb"
          },
          tooltip: {
            pointFormat: '{point.x}'
          },
          data: data_scatter
        }
      ],
      credits: {
        enabled: false
      },
      exporting: {
        enabled: false
      }
 },
 function(chart) {
  chart.xAxis[0].plotLinesAndBands[0].label.toFront();
});
}
/* RIDGLINE DENSITY PLOT */

function processDensity(step, precision, densityWidth, data) {
  let xiData = [];
  
  function prcessXi(args) {
    let tempXdata = [];
    let tileSteps = 3; // Nbr of point at the top and end of the density
    let min = Infinity, max = -Infinity;

      //process the range of the data set
      min = Math.min(min, Math.min(...args));
      max = Math.max(max, Math.max(...args));

    for (i = min - tileSteps * step; i < max + tileSteps * step; i++) {
      tempXdata.push(i);
    }
    return tempXdata;
  }

  data.forEach((d) => {
    xiData = prcessXi(d);
  })

  let gap = -1;

  function density(dataSource) {
    gaussian_bw = 0.35;
    
    let data = [];
    let N = dataSource.length;

    gap++;
    for (i = 0; i < xiData.length; i++) {
      let temp = 0;
      for (j = 0; j < dataSource.length; j++) {
        temp = temp + GaussKDE(xiData[i], dataSource[j], gaussian_bw);
      }
      data.push([xiData[i], (1 / N) * temp]);
    }

    return data.map((densityPoint, i) => {

      if (densityPoint[1] > precision) {
        return [xiData[i], gap, densityPoint[1] * densityWidth + gap];
      } else {
        return [xiData[i], null, null];
      }
    });
  }

  let results = [];
  let index = 0;

  data.forEach((e) => {
      results.push([]);
      results[index] = density(e).slice();
      index++;
    });

    return { xiData, results };
}

function createRidgelineDensityPlot(background_csv, tinyt_csv, thresholds_csv) {

  var data_background = processCSV(background_csv);
  var data_sample = processCSV(tinyt_csv);
  
  var obj_background = []
  
  data_background.forEach(function(d) {
    obj_background.push({Deletion: d.Deletion, ScaledProportion: d.ScaledProportion})
  });
  

	var	obj_sample = []

	data_sample.forEach(function(d) {
		obj_sample.push({Deletion: d.Deletion, ScaledProportion: d.ScaledProportion})
	});

  
  let deletions_list = Object.values(unique_deletions).reverse() // Take the deletions list and reverse so the smallest value would be on top
    
    // Create data arrays from BACKGROUND object values and update the data
    let dataArray = [];
    for (i = 0; i < deletions_list.length; i++) {
      dataArray.push([]);
    }
    obj_background.forEach((e) => {
      deletions_list.forEach((key, value) => {
        if (e.Deletion == key) {
          dataArray[value].push(e.ScaledProportion);
        }
      });
    });

    let dataArray_zscore = []
    for (i = 0; i < deletions_list.length; i++) {
      dataArray_zscore.push([]);
    }

    dataArray.forEach((d, i) => {
      let backgr_mean = MathAverage(d)
      let backgr_sd = MathStandardDeviation(d)
      let Zscore = Infinity

      d.forEach(function(e) {
        Zscore = MathZscore(e, backgr_mean, backgr_sd)
        dataArray_zscore[i].push(Zscore)
      })
    });
  

  // Create data arrays from SAMPLE object values and update the data
  let sampleDataArray = [];
    for (i = 0; i < deletions_list.length; i++) {
      sampleDataArray.push([]);
    }
    obj_sample.forEach((e) => {
      deletions_list.forEach((key, value) => {
        if (e.Deletion == key) {
          sampleDataArray[value].push(e.ScaledProportion);
        }
      });
    });

    let sampleDataArray_zscore = []
    for (i = 0; i < deletions_list.length; i++) {
      sampleDataArray_zscore.push([]);
    }
    dataArray.forEach((d, i) => { // Use background data to calculate the mean like previously
      let backgr_mean = MathAverage(d)
      let backgr_sd = MathStandardDeviation(d)
      let Zscore = Infinity

      Zscore = MathZscore(sampleDataArray[i], backgr_mean, backgr_sd)
      sampleDataArray_zscore[i].push(Zscore)
    });

    // Process density data
    let step = 1.5, width = 1.3

    let data = processDensity(step, precision, width, dataArray_zscore)

    // Structure the data to create the chart
    let chartsNbr = data.results.length
    let xi = data.xiData

    // Add BACKGROUND data to series
    let dataSeries = [], series = [];
    data.results.forEach((e, i) => {
      dataSeries.push([]);
      dataSeries[i] = e;
      series.push({
          name: 'Background',
          tooltip: {
            pointFormat: 'Deletion: ' + deletions_list[i] + '<br />Z-score: {point.x}'
          },
          connectNulls: true,
          zIndex: chartsNbr - i,
          data: dataSeries[i],
          color: Highcharts.getOptions().colors[i%10] // Take only the last number since we have only 10 color values specified
      });
    });

  // Add SAMPLES data to series
  sampleDataArray_zscore.forEach((zscore, i) => {
    series.push({
          type: 'scatter',
          name: 'Analysed sample',
          marker: {
            enabled: true,
            symbol: "diamond",
            fillColor: Highcharts.getOptions().colors[i%10], // Take only the last number since we have only 10 color values specified
            height: "10px",
            lineWidth: (zscore >= 3 ? 4 : '0'),
						lineColor: (zscore >= 3 ? "#F10E43" : '')

        },
				tooltip: {
					pointFormat: 'Deletion: ' + deletions_list[i] + '<br />Z-score: {point.x}'
				},
        zIndex: chartsNbr + i,
        data: [{
            x: zscore, // X axis value is the zscore
            y: i + 0.2, // Move the marker bit higher from bottom line
        }]
    })
  });

    Highcharts.chart("plot_ridgeline_density", {
      chart: {
        type: "areasplinerange",
        height: "460px", // Make it higher
        animation: true,

      },
      name: 'Background',
      title: {
        text: "All deletions for " + gene_id,
        style: { "fontSize": "14px" }


      },
      xAxis: {
        title: { text: "Z-score" },
        labels: { format: "{value}" },
        tickInterval: 1
      },
      yAxis: {
        title: { text: null },
        categories: deletions_list,
        max: data.results.length,
        labels: {
          formatter: function () {
            if (this.pos < chartsNbr) return this.value;
          },
          style: {
            textTransform: "capitalize",
            fontSize: "9px"
          }
        },
        startOnTick: true,
        gridLineWidth: 1,
        tickmarkPlacement: "on"
      },
      plotOptions: {
        areasplinerange: {
          marker: {
            enabled: false
          },
          states: {
            hover: {
              enabled: false
            }
          },
          pointStart: xi[0],
          animation: {
            duration: 0
          },
          fillOpacity: 0.1,
          lineWidth: 1
        }
      },
      legend: {
        enabled: false
      },
      series: series,
      credits: {
        enabled: false
      },
      exporting: {
        enabled: false
      }
    });

    return sampleDataArray_zscore // Function returns array of Z-scores
}
  