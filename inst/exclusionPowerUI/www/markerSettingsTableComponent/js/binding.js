let binding = new Shiny.InputBinding();

const buildMutationModelSelect = function(marker) {
  const mutModels = ['None'];
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

const buildRow = function(marker) {
  let sexLinked;
  if (marker.chrom == '23')
    sexLinked = '1';
  else
    sexLinked = '0';

  return '<tr data-marker="' + marker.markerName + '">' +
    '<td>' + marker.markerName + '</td>' +
    '<td>' + buildCheckbox('sexLinked-checkbox', marker.chrom == '23') + '</td>' +
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
      let chrom = $(row).find('input.sexLinked-checkbox')[0].checked;

      // translate to pedtools chromosome index
      if (chrom)
        chrom = 23;
      else
        chrom = 'NA';

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
