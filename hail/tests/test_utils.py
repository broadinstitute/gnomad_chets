import unittest

from utils import *

hc = None
verbose = False


def setUpModule():
    global hc
    global verbose
    verbose = '-v' in sys.argv
    hc = HailContext()  # master = 'local[2]')


def tearDownModule():
    global hc
    hc.stop()
    hc = None


def test_vds_from_rows(rows, flat_schema, types):
    """

    :param rows: Data
    :param list of str flat_schema: Names for schema
    :param list of obj types: Types for schema
    :return: Fake VDS
    :rtype: VariantDataset
    """
    flat_schema.insert(0, 'v')
    types.insert(0, TVariant())
    schema = TStruct(flat_schema, types)
    for i, row in enumerate(rows):
        row['v'] = Variant.parse('1:{}:A:T'.format(i + 1000))
    return VariantDataset.from_table(KeyTable.parallelize(rows, schema, key='v'))


class FilteringTests(unittest.TestCase):

    @staticmethod
    def create_filter_test_vds():
        """

        :return: VDS with some filters
        :rtype: VariantDataset
        """
        rows = [
            # Bi-allelic expected behavior
            {'v': Variant.parse('1:10000:A:T'),   'InbreedingCoeff': None, 'AS_FilterStatus': [[]],              'expected_filters': [],                        'expected_after_split': [[]]},
            {'v': Variant.parse('1:10001:A:T'),   'InbreedingCoeff': 0.0,  'AS_FilterStatus': [[]],              'expected_filters': [],                        'expected_after_split': [[]]},
            {'v': Variant.parse('1:10002:A:T'),   'InbreedingCoeff': -0.5, 'AS_FilterStatus': [[]],              'expected_filters': ['InbreedingCoeff'],       'expected_after_split': [['InbreedingCoeff']]},
            {'v': Variant.parse('1:10003:A:T'),   'InbreedingCoeff': -0.5, 'AS_FilterStatus': [['RF']],          'expected_filters': ['InbreedingCoeff', 'RF'], 'expected_after_split': [['InbreedingCoeff', 'RF']]},
            {'v': Variant.parse('1:10004:A:T'),   'InbreedingCoeff': 0.0,  'AS_FilterStatus': [['RF']],          'expected_filters': ['RF'],                    'expected_after_split': [['RF']]},
            {'v': Variant.parse('1:10005:A:T'),   'InbreedingCoeff': 0.0,  'AS_FilterStatus': [['RF', 'AC0']],   'expected_filters': ['RF', 'AC0'],             'expected_after_split': [['RF', 'AC0']]},

            # Multi-allelic expected behavior
            {'v': Variant.parse('2:10000:A:T,C'), 'InbreedingCoeff': 0.0,  'AS_FilterStatus': [[], []],          'expected_filters': [],                               'expected_after_split': [[], []]},
            {'v': Variant.parse('2:10001:A:T,C'), 'InbreedingCoeff': -0.5, 'AS_FilterStatus': [[], []],          'expected_filters': ['InbreedingCoeff'],              'expected_after_split': [['InbreedingCoeff'], ['InbreedingCoeff']]},
            {'v': Variant.parse('2:10002:A:T,C'), 'InbreedingCoeff': 0.0,  'AS_FilterStatus': [['RF'], []],      'expected_filters': [],                               'expected_after_split': [['RF'], []]},
            {'v': Variant.parse('2:10003:A:T,C'), 'InbreedingCoeff': 0.0,  'AS_FilterStatus': [['RF'], ['RF']],  'expected_filters': ['RF'],                           'expected_after_split': [['RF'], ['RF']]},
            {'v': Variant.parse('2:10004:A:T,C'), 'InbreedingCoeff': 0.0,  'AS_FilterStatus': [['RF'], ['AC0']], 'expected_filters': ['RF', 'AC0'],                    'expected_after_split': [['RF'], ['AC0']]},
            {'v': Variant.parse('2:10005:A:T,C'), 'InbreedingCoeff': -0.5, 'AS_FilterStatus': [['RF'], []],      'expected_filters': ['InbreedingCoeff'],              'expected_after_split': [['InbreedingCoeff', 'RF'], ['InbreedingCoeff']]},
            {'v': Variant.parse('2:10006:A:T,C'), 'InbreedingCoeff': -0.5, 'AS_FilterStatus': [['RF'], ['AC0']], 'expected_filters': ['InbreedingCoeff', 'RF', 'AC0'], 'expected_after_split': [['InbreedingCoeff', 'RF'], ['InbreedingCoeff', 'AC0']]},

            # Unexpected behavior
            {'v': Variant.parse('9:10000:A:T'),   'InbreedingCoeff': 0.0,  'AS_FilterStatus': None,              'expected_filters': None,                      'expected_after_split': None},
            {'v': Variant.parse('9:10001:A:T'),   'InbreedingCoeff': None, 'AS_FilterStatus': None,              'expected_filters': None,                      'expected_after_split': None},
            {'v': Variant.parse('9:10002:A:T'),   'InbreedingCoeff': -0.5, 'AS_FilterStatus': None,              'expected_filters': None,                      'expected_after_split': None},
            {'v': Variant.parse('9:10003:A:T'),   'InbreedingCoeff': 0.0,  'AS_FilterStatus': [None],            'expected_filters': None,                      'expected_after_split': None},
            {'v': Variant.parse('9:10004:A:T'),   'InbreedingCoeff': 0.0,  'AS_FilterStatus': [[None]],          'expected_filters': None,                      'expected_after_split': None},
            {'v': Variant.parse('9:10005:A:T,C'), 'InbreedingCoeff': 0.0,  'AS_FilterStatus': [[], None],        'expected_filters': None,                      'expected_after_split': [[], None]},
            {'v': Variant.parse('9:10006:A:T,C'), 'InbreedingCoeff': 0.0,  'AS_FilterStatus': [[], [None]],      'expected_filters': None,                      'expected_after_split': [[], None]},
            {'v': Variant.parse('9:10007:A:T,C'), 'InbreedingCoeff': 0.0,  'AS_FilterStatus': [['RF'], [None]],  'expected_filters': None,                      'expected_after_split': [['RF'], None]},
        ]
        schema = ['v', 'InbreedingCoeff', 'AS_FilterStatus', 'expected_filters', 'expected_after_split']
        types = [TVariant(), TDouble(), TArray(TSet(TString())), TSet(TString()), TArray(TSet(TString()))]
        return VariantDataset.from_table(KeyTable.from_py(hc, rows, TStruct(schema, types), key_names=['v']))

    @classmethod
    def setUpClass(cls):
        cls.vds = cls.create_filter_test_vds()
        if verbose: cls.vds.variants_table().show(50)

    def test_allele_filtering(self):
        site_filters = {
            'InbreedingCoeff': 'isDefined(va.InbreedingCoeff) && va.InbreedingCoeff < -0.3'
        }

        result_vds = set_site_filters(self.vds, site_filters, 'va.AS_FilterStatus')
        if verbose: result_vds.variants_table().show(50)
        result = result_vds.query_variants('variants.map(v => (isMissing(va.filters) && isMissing(va.expected_filters)) || va.filters == va.expected_filters).counter()')
        self.assertEqual(result[True], sum(result.values()))

        split_vds = result_vds.split_multi().annotate_variants_expr(index_into_arrays(['va.AS_FilterStatus', 'va.expected_after_split']))
        result_split_vds = set_site_filters(split_vds, site_filters, 'va.AS_FilterStatus')
        if verbose: result_split_vds.variants_table().show(50)
        result = result_split_vds.query_variants('variants.map(v => (isMissing(va.filters) && isMissing(va.expected_filters)) || va.filters == va.expected_after_split).counter()')
        self.assertEqual(result[True], sum(result.values()))


class KeyTableTests(unittest.TestCase):

    @staticmethod
    def create_frequency_kt():
        """
        KeyTable with some frequency data

        :return: keytable with frequency data
        :rtype: KeyTable
        """
        rows = [
            # Bi-allelic expected behavior
            {'v': Variant.parse('1:10000:A:T'), 'AC_NFE': 1,  'AC_AFR': 8,   'Hom_NFE': 0, 'Hom_AFR': 0},
            {'v': Variant.parse('1:10001:A:T'), 'AC_NFE': 10, 'AC_AFR': 100, 'Hom_NFE': 1, 'Hom_AFR': 10},
        ]
        schema = ['v', 'AC_NFE', 'AC_AFR', 'Hom_NFE', 'Hom_AFR']
        types = [TVariant(), TInt(), TInt(), TInt(), TInt()]
        return KeyTable.from_py(hc, rows, TStruct(schema, types), key_names=['v'])

    @classmethod
    def setUpClass(cls):
        cls.kt = cls.create_frequency_kt()
        if verbose: cls.kt.show(50)

    def test_melt_kt(self):
        melted_kt = melt_kt(self.kt, columns_to_melt=['AC_NFE', 'AC_AFR', 'Hom_NFE', 'Hom_AFR'])
        self.assertEqual(melted_kt.count(), 8)
        self.assertEqual(sorted(melted_kt.columns), sorted(['v', 'value', 'variable']))
        self.assertEqual(melted_kt.query('variable.counter()'), {'AC_NFE': 2, 'AC_AFR': 2, 'Hom_NFE': 2, 'Hom_AFR': 2})
        self.assertEqual(melted_kt.filter('v == Variant("1:10000:A:T") && variable == "AC_NFE"').query('value.collect()'), [1])
        if verbose: melted_kt.show(50)

    def test_melt_grouped_kt(self):
        grouped_melted_kt = melt_kt_grouped(self.kt,
                                            columns_to_melt={'NFE': ['AC_NFE', 'Hom_NFE'], 'AFR': ['AC_AFR', 'Hom_AFR']},
                                            value_column_names=['AC', 'Hom'],
                                            key_column_name='pop')
        self.assertEqual(grouped_melted_kt.count(), 4)
        self.assertEqual(sorted(grouped_melted_kt.columns), sorted(['v', 'pop', 'AC', 'Hom']))
        self.assertEqual(grouped_melted_kt.query('pop.counter()'), {'NFE': 2, 'AFR': 2})
        if verbose: grouped_melted_kt.show(50)



if __name__ == '__main__':
    unittest.main()