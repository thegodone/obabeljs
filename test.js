'use strict';

var ob = require('./dist/obabel.js');

console.log(ob.hello);

var mol = ob.Mol.fromSmiles('CCCCC');


var canindex = mol.canonicalindex();


console.log(canindex.size())

mol.delete();

/*describe('OBabel loading', function () {
    it('JS functions Hello', function () {
        ob.hello.should.equal('world');

    });
});
*/