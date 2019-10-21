const buildDropdown = function(field, selected) {
  return '<select name="' + field.name + '">' +
  field.options.map(function(opt) {
    let selectedAttr = (opt.value == selected) ? ' selected' : '';
    return '<option value="' + opt.value + '"' + selectedAttr + '>' + opt.label + '</option>';
  }).join(' ') +
  '</select>';
};

const buildCheckbox = function(field, checked) {
  let selectedAttr = checked ? ' checked="checked"' : '';
  return '<input type="checkbox" name="' + field.name + '"' + selectedAttr + '>';
};

const buildRadioButtons = function(field, selected) {
  return '<form class="form-inline">' +
  field.options.map(function(opt) {
    let selectedAttr = (opt.value == selected) ? ' checked="checked"' : '';
    return '<label class="radio-inline">' +
           '<input type="radio" ' + selectedAttr +
                   'name="' + field.name + '" ' +
                   'value="' + opt.value + '">' +
                   opt.label + '</label>';
  }).join(' ') + '</form>';
};

const buildHeader = function(fields) {
  return '<thead>' +
  fields.map(function(field) {
    return '<th>' + field.name + '</th>';
  }).join(' ') +
  '</thead>';
};

const buildRow = function(row, fields) {
  return '<tr>' +
  fields.map(function(field) {
    let ret = '<td data-input-type="' + field.type + '" data-input-name="' + field.name + '">';
    if (field.type == 'checkbox') {
      ret += buildCheckbox(field, row[field.name][0]);
    } else if (field.type == 'radio') {
      ret += buildRadioButtons(field, row[field.name]);
    } else if (field.type == 'dropdown') {
      ret += buildDropdown(field, row[field.name]);
    } else {
      ret += row[field.name];
    }
    return ret + '</td>';
  }).join(' ') + '</tr>';
};

const buildTable = function(el, state, fields) {
    let table = $(el).find('table');
    table.empty();
    table.append(buildHeader(fields));

    table.append('<tbody></tbody>');
    let tbody = table.find('tbody');
    state.forEach(function(row) {
      tbody.append(buildRow(row, fields));
    });
};

const parseTable = function(el) {
  return Array.prototype.map.call($(el).find('tbody tr'), function(row) {
    let record = {};
    $(row).find('td').each(function(i, field) {
      let type = $(field).data('input-type');
      let name = $(field).data('input-name');
      switch(type) {
        case 'checkbox':
          record[name] = $(field).find('input[type="checkbox"]')[0].checked;
          break;

        case 'radio':
          record[name] = $(field).find('input:radio:checked').val();
          break;

        case 'dropdown':
          record[name] = $(field).find('select').val();
          break;

        default:
          record[name] = field.innerText;
      }
    });

    return record;
  });
};

let binding = new Shiny.InputBinding();

$.extend(binding, {
  find: function(scope) {
    return $(scope).find('.settings-table');
  },

  initialize: function(el) {
  },

  subscribe: function(el, callback) {
    $(el).on('change', 'input, select', callback);
  },

  receiveMessage: function(el, message) {
    buildTable(el, message.data, message.fields);
  },

  getValue: function(el) {
    return JSON.stringify(parseTable($(el).find('table')));
  }
});

Shiny.inputBindings.register(binding);
