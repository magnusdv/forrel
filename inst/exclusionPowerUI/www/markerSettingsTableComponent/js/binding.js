let binding = new Shiny.InputBinding();

const mutModels = ['None'];

const buildMutationModelSelect = function(marker) {
  let ret = '<select class="mutationModel-select">';
  mutModels.forEach(function(model) {
    if (model == marker.mutationModel)
      ret += '<option selected value="' + model + '">' + model + '</option>';
    else
      ret += '<option value="' + model + '">' + model + '</option>';
  });
  ret += '</select>';
  return ret;
};

const buildCheckbox = function(klass, selected) {
  if (selected)
    return '<input checked type="checkbox" class="' + klass + '">';
  else
    return '<input type="checkbox" class="' + klass + '">';
};

const buildRadioButtons = function(klass, options, selected) {

  return options.map(function(opt) {
    selectedAttr = '';
    if (selected == opt.value) selectAttr = ' selected ';
    return '<label class="radio-inline">' +
           '<input type="radio" ' + selectedAttr +
                   'name="' + opt.name + '" ' +
                   'value="' + opt.value + '" ' +
                   'class="' + klass + '">' +
                   opt.label + '</label>';
  }).join(' ');
};

const sexLinkedRadioOptionsForMarker = function(markerName) {
  return [
    {
      name: markerName + '-radio',
      value: '1',
      label: 'Autosomal'
    },
    {
      name: markerName + '-radio',
      value: '23',
      label: 'X'
    }
  ];
};

const buildRow = function(marker) {
  let sexLinkedSelected = (marker.chrom == '23') ? '23' : '1';
  return '<tr data-marker="' + marker.markerName + '">' +
    '<td>' + marker.markerName + '</td>' +
    '<td>' + buildRadioButtons('sexLinked-radio',
                               sexLinkedRadioOptionsForMarker(marker.markerName),
                               sexLinkedSelected) + '</td>' +
    '<td>' + buildMutationModelSelect(marker) + '</td>' +
    '<td>' + buildCheckbox('includeInCalculation-checkbox', marker.includeInCalculation == '1') + '</td>' +
    '</tr>';
};

const buildTable = function(el, state) {
    $(el).find('tbody').empty();

    state.forEach(function(marker) {
      $(el).find('tbody').append(buildRow(marker));
    });
};

$.extend(binding, {
  find: function(scope) {
    return $(scope).find('.marker-settings-table');
  },

  initialize: function(el) {
    console.log('initialize called for marker-settings-table widget');
  },

  subscribe: function(el, callback) {
    console.log('subscribe called for marker-settings-table widget');
    $(el).on('change', 'input, select', callback);
  },

  receiveMessage: function(el, message) {
    console.log('messageReceived');
    console.log(message);

    buildTable(el, message);
  },

  getValue: function(el) {
    console.log('trying to get widget value');

    let res = [];

    $(el).find('tbody tr').each(function callback(_, row) {
      let markerName = $(row).data('marker');
      let includeInCalculation = $(row).find('input.includeInCalculation-checkbox')[0].checked;
      let chrom = $(row).find('.sexLinked-radio:checked').val();

      let mutationModel = $(row).find('select.mutationModel-select').val();

      res.push({
        markerName: markerName,
        includeInCalculation: includeInCalculation,
        chrom: chrom,
        mutationModel: mutationModel
      });
    });

    console.log(res);

    return JSON.stringify(res);
  },
});

Shiny.inputBindings.register(binding);
