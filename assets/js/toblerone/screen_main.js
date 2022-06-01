/* THIS SECTION CONTAINS CODE FOR THE INTERFACE FUNCTIONALITY */
var gene_id;

window.addEventListener('pywebviewready', function() {
    setTimeout(function() {
		pywebview.api.getListOfGenesAndCohorts().then(updateMainSelection) // Get list of genes and cohorts from the catalogue
	}, 100);

	setTimeout(function() {
		pywebview.api.getDynamicConfiguration().then(fillConfigurationValues) // Update configuration section
	}, 200);
})

$('.dropdown-toggle').dropdown()

var listOfGenesCohorts;

// Update genes selection on main page
function updateMainSelection(res) {
	listOfGenesCohorts = res

	$.each(res.Genes, function(i, element) {
		$("#gene").append($('<option>', { 
			value: element.Gene,
			text : element.Gene 
		}))
	})
}

function fillConfigurationValues(conf) {
	$("#config_DeleteTempFiles").prop("checked", conf["DeleteTempFiles"]);
	$("#config_TwoPass").prop("checked", conf["TwoPass"]);
}

 
// Validate form
window.addEventListener('load', () => {
  var forms = document.getElementsByClassName('needs-validation');
  for (let form of forms) {
    form.addEventListener('submit', (evt) => {
      if (!form.checkValidity()) {
        evt.preventDefault();
        evt.stopPropagation();
      } else {
        evt.preventDefault();
		startAnalysis() // Run the analysis if form was validated
      }
      form.classList.add('was-validated');
    });
  }
})


// Open file open dialog box
function selectFile() {
	pywebview.api.selectSeqFile().then(setFilePath)
}

// When BAM file is selected
function setFilePath(filePath) {
	document.getElementById("bamfile").value = filePath
	if (filePath) {
		$("#bamfilelabel").html(filePath)
		if ($('#refgenome').val() != "") {
			changeInputIcon("refgenome")
		}
	}
	else {
		$("#bamfilelabel").html('Choose a BAM file')
	}
}

// Change icons (e.g. when selections are made then display a checkmark icon)
function changeInputIcon(inputId, type) {
	if (type == 'check') {
		$('#' + inputId).parent().parent().find('.stepper-icon').removeClass('bg-secondary')
		$('#' + inputId).parent().parent().find('.stepper-icon').addClass('bg-primary')
		$('#' + inputId).parent().parent().find('.stepper-icon').children('.material-icons').html("check")
	}
	
	if (type == 'edit') {
		$('#' + inputId).parent().parent().find('.stepper-icon').removeClass('bg-primary')
		$('#' + inputId).parent().parent().find('.stepper-icon').addClass('bg-secondary')
		$('#' + inputId).parent().parent().find('.stepper-icon').children('.material-icons').html("edit")
	}
}

// Iniate the analysis at backend
function startAnalysis() {
	var cohort = $('[name=cohort]').val()
	var gene = $('[name=gene]').val()
	var refgenome = $('[name=refgenome]').val()
	var bamfile = $('[name=bamfile]').val()

	gene_id = gene;

	if (cohort && gene && refgenome && bamfile) {
		$("#progressbar").removeClass('d-none');
		pywebview.api.analyseSample(gene, cohort, refgenome, bamfile).then(showResponse)
	}
}


// If gene is selected or changed, then update the cohorts selection
$("#gene").change(function() {
	if ($(this).data('options') === undefined) {
	  $(this).data('options', $('#cohort option').clone());
	}

	gene_id = $(this).val()

	$("#cohort").html('')
	$.each(listOfGenesCohorts.Cohorts, function(i, element) {
		if (element.Gene == gene_id) {
			$("#cohort").append($('<option>', { 
				value: element.ID,
				text : element.Name 
			}));
		}
	})
  })

// If settings are saved then send them to backend and close the settings section
$("#btn_save_settings").click(function() {
	let del_temp_files_val = $("#config_DeleteTempFiles").is(":checked")
	let twopass_val = $("#config_TwoPass").is(":checked")

	let new_conf = {
		"DeleteTempFiles": del_temp_files_val,
		"TwoPass": twopass_val
	}
	
	pywebview.api.updateConfiguration(new_conf).then(function() {
		$("#toggleSettings").collapse('toggle') 
	})
})

// If option in select box is changed
$("select").on('change', function(e) {
	if (e.target.id == "refgenome" && $('#bamfile').val() == "") {
		return;
	}
		changeInputIcon(e.target.id, "check")
})

// If a (different) Gene is selected
$("#gene").on('change', function(e) {
	if ($('#gene').val()) {
		changeInputIcon("gene", "check")
		changeInputIcon("cohort", "check")
	}
	else {
		changeInputIcon("gene", "edit")
		changeInputIcon("cohort", "edit")
	}
})
