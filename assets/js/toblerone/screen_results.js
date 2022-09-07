/* THIS SECTION CONTAINS CODE WHICH IS USED FOR THE RESULTS MODAL TO FUNCTION, PLOTS ARE EXCLUDED FROM HERE */


var analysis_info;
var tinyt_csv;
var background_csv;
var thresholds_csv;
var deletion_id;
var selected_data_type = "proportion"; // Default data type is porportion (instead of zscore)
var tab1 = 'res-summary-content'
var tab2 = 'res-density-content'
var unique_deletions = {};

// Fill textarea with error text and display it to user
function errorHandling(response) {
	$("#raw-stderr").val(response.tinyt_results)

	$("#progressbar").addClass('d-none')
	$('#resultsModal').modal('show')
	$('.nav-tabs a[href="#res-stderr-content"]').tab('show');
}

// When receiving a response from backend then either show the results or an error
function showResponse(response) {
	// If tinyt returns an error, then handle that and break the process
	if (response.errors.tinyt) {
		errorHandling(response)
		return;
	}

	// Set global variable contents
	tinyt_csv = response.tinyt_results
	tinyt_log = response.tinyt_log
	background_csv = response.background
	thresholds_csv = response.thresholds
	analysis_info = response.analysis_info

	first_deletion_id = updateDeletionSelection(background_csv)
	deletion_id = first_deletion_id // Set the first deletion ID which will be then used to show the first plot

	$("#raw-tinyt").val(tinyt_csv)
	$("#raw-stderr").val(tinyt_log)

	showDataOnSummaryTab() // Show data on the first summary tab
	filterDataByDeletion(deletion_id, selected_data_type) // For other data/plots, filter data by deletion name
}
	
// Update data under the Deletion tab
function updateDeletionSelection(data) {
	var data_arr = csvToArray(data);

	$.each(data_arr, function(i, item) {
		if (item.Deletion) {
			unique_deletions[item.Deletion] = item.Deletion;
		}
	});

	$("#deletion_selection").empty() // Empty it first to avoid duplicates

	$.each(unique_deletions, function(i, element) {
		$("#deletion_selection").append($('<option>', { 
			value: element,
			text : element
		}))
	})

	first_element = unique_deletions[Object.keys(unique_deletions)[0]]

	$("#deletion_selection").val(first_element).select() // Select the first element on the selection menu

	return first_element // Return first deletion ID
}


// Use deletion name to filter the data
function filterDataByDeletion(deletion, data_type) {
	deletion_id = deletion;
	var data_thresholds = processThresholdsCSV(thresholds_csv, deletion) // Get thresholds for this deletion_id 

	// Give data a values from the csv, separate columns by spaces
	var data_background = processCSV(background_csv, deletion);
	var data_sample = processCSV(tinyt_csv, deletion)[0];
	var sample_value = { 
		ScaledProportion: data_sample.ScaledProportion,
		Zscore: null
	}
	var sample_threshold_data = getSamplesThresholdData(data_thresholds, sample_value.ScaledProportion)


	// If zscore is selected from data type to display, then take background ScaledProportions into an array and calculate z score from samples ScaledProportion and background data
	if (data_type == "zscore") {
		var backgr_scaledprop_arr = []
		data_background.forEach(function(d) {
			backgr_scaledprop_arr.push(d.ScaledProportion)
		});

		backgr_mean = MathAverage(backgr_scaledprop_arr)
		backgr_sd = MathStandardDeviation(backgr_scaledprop_arr)

		// Update background data
		data_background.forEach(function(d) {
			d.Zscore = MathZscore(d.ScaledProportion, backgr_mean, backgr_sd)
		});

		// Add zscore to samples value, calculated from ScaledProportion
		sample_value["Zscore"] = MathZscore(sample_value.ScaledProportion, backgr_mean, backgr_sd)
	}

	createDensityPlot(data_background, sample_value, data_thresholds, data_type, deletion) // Create the density plot

    showDataOnDensityTab(data_sample, sample_value, sample_threshold_data, data_background, data_type) // Show text data on tab

	$("#progressbar").addClass('d-none')
	$('#resultsModal').modal('show')
}

// Show data on the first summary tab
function showDataOnSummaryTab() {
	// Set the boxplot fixed to del4_5_6_7 instead of being dynamic
	let fix_deletion_id = "del4_5_6_7"
	let data_background_del4_5_6_7 = processCSV(background_csv, fix_deletion_id)
	let data_sample_del4_5_6_7 = processCSV(tinyt_csv, fix_deletion_id)[0];
	let data_thresholds = processThresholdsCSV(thresholds_csv, fix_deletion_id)
	let sample_value = data_sample_del4_5_6_7.ScaledProportion
    let sample_threshold_data_del4_5_6_7 = getSamplesThresholdData(data_thresholds, sample_value) // del4-7 detection status
	createBoxplot(data_background_del4_5_6_7, data_sample_del4_5_6_7, fix_deletion_id)

	// Display the ridgeline density plot
	var zscores_all_dels = createRidgelineDensityPlot(background_csv, tinyt_csv, thresholds_csv) // Returns Z-scores for all deletions

	let del_data_for_summary = sample_threshold_data_del4_5_6_7 // Assign the del4-7 data to show in summary tab.
	let outliers_total = zscores_all_dels.filter(x => x >= 3).length // Count how many has Zscore >= 3
	let summary_text_w_values_text = `The ${del_data_for_summary.Gene} ${del_data_for_summary.Deletion} was ${del_data_for_summary.Label} (shown on left). There were ${outliers_total} total outliers in IKZF1 deletions (on the right). See Deletions for more details.` // Create the text
	$('#summary_text_w_values').html(summary_text_w_values_text) // Push it to the HTML

	let thresholdLabel = sample_threshold_data_del4_5_6_7.Label
	let textColour = sample_threshold_data_del4_5_6_7.Colour
	let metadata = analysis_info

	$('#' + tab1 + ' .info_sample').html(metadata.SampleFile)
	$('#' + tab1 + ' .info_gene').html(metadata.Gene)
	$('#' + tab1 + ' .info_coordinates').html(metadata.GeneCoordinates)
	$('#' + tab1 + ' .info_genelength').html(data_sample_del4_5_6_7.GeneLength)
	$('#' + tab1 + ' .info_cohort').html(metadata.CohortName)
	$('#' + tab1 + ' .results_clone').html(thresholdLabel)
	$('#' + tab1 + ' .results_clone').css('color', textColour)
	$('#' + tab1 + ' .results_clone').css('border', '1px solid ' + textColour)
}

// Sow data on the density tab
function showDataOnDensityTab(sample_data, sample_value, sample_threshold, data_background, data_type) {
	let thresholdLabel = sample_threshold.Label
	let textColour = sample_threshold.Colour

	// If zscore is selected from data type to display, then take background ScaledProportions into an array and calculate z score from samples ScaledProportion and background data
	if (data_type == "zscore") {
		$('#' + tab2 + ' .results_scaledproportions_title').html("Z score")
		$('#' + tab2 + ' .results_scaledproportions').html(sample_value.Zscore)
	}
	else {
		$('#' + tab2 + ' .results_scaledproportions_title').html("Scaled Proportion")
		$('#' + tab2 + ' .results_scaledproportions').html(sample_value.ScaledProportion)
	}

	$('#' + tab2 + ' .results_clone').html(thresholdLabel)
	$('#' + tab2 + ' .results_gene').html(sample_data.Gene)
	$('#' + tab2 + ' .results_deletion').html(sample_data.Deletion)
	$('#' + tab2 + ' .results_count').html(sample_data.Count)
	$('#' + tab2 + ' .results_total').html(sample_data.Total)
	$('#' + tab2 + ' .results_genelength').html(sample_data.GeneLength)
	$('#' + tab2 + ' .results_readlength').html(sample_data.ReadLength)
	$('#' + tab2 + ' .results_clone').css('color', textColour)
	$('#' + tab2 + ' .results_clone').css('border', '1px solid ' + textColour)
}

/* THIS SECTION CONTAINS CODE FOR PROCESSING DATA BEFORE MAKING PLOTS */

// Convert CSV to array and filter data by deletion ID
function processThresholdsCSV(data, deletion_id) {
    let data_arr = csvToArray(data);

    data_arr_filtered = data_arr.filter(function(item) {
		return item.Deletion == deletion_id
	})

	if (data_arr_filtered.length == 0) {
		data_arr_filtered = data_arr.filter(function(item) {
			return item.Deletion == '*'
		})
	}

    return data_arr_filtered;
}

// Process both background and sample CSV file -> convert to array and filter out all other besides the one given deletion ID
function processCSV(data, deletion_id = null) {
    let data_arr = csvToArray(data);

	if (deletion_id) {
		data_arr = data_arr.filter(function(item) {
			return item.Deletion == deletion_id
		})
	}

    data_arr.forEach(function(d) {
		d.Count = +d.Count
		d.ScaledProportion = +d.ScaledProportion
	});

    return data_arr;
}

// Convert CSV to array
function csvToArray(str, delimiter = ",") {
	// Converts a CSV file/content into an array of objects
	const headers = str.slice(0, str.indexOf("\n")).split(delimiter);

	const rows = str.slice(str.indexOf("\n") + 1).split("\n");

	const arr = rows.map(function (row) {
	  const values = row.split(delimiter);
	  const el = headers.reduce(function (object, header, index) {
		object[header] = values[index];
		return object;
	  }, {});
	  return el;
	});

	return arr;
}

/* Mathematical functions such as average, SD, Zscore and also for creating range of numbers */
const MathAverage = (array) => array.reduce((a, b) => a + b) / array.length;

const MathStandardDeviation = (arr, usePopulation = true) => {
	let arr_mean = MathAverage(arr)
	return Math.sqrt(
	  arr
		.reduce((acc, val) => acc.concat((val - arr_mean) ** 2), [])
		.reduce((acc, val) => acc + val, 0) /
		(arr.length - (usePopulation ? 0 : 1))
	);
  };

  const MathZscore = (x, mean, sd) => {
	  if ((x - mean) == 0) { // Return 0 if we cannot calculate
		  return 0
	  }
	  else {
		return (x - mean) / sd
	}
  };

  const MathRange = (min, max) => [...Array(max - min + 1).keys()].map(i => i + min);


// Capture events when user changes deletion on select menu
$("#deletion_selection").on('change', function(e) {
	filterDataByDeletion(e.target.value, selected_data_type)
})

// Capture events when user changes data type on select menu
$("#data_selection").on('change', function(e) {
	selected_data_type = e.target.value
	filterDataByDeletion(deletion_id, e.target.value)
})