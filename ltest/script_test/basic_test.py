# Test script for genome_util package - it should be launched from
# the root of the genome_util module, ideally just with 'make test', as
# it looks for a hardcoded relative path to find the 'test.cfg' file
import unittest
import json
import ConfigParser

from pprint import pprint

from subprocess import call

from biokbase.auth import Token

# Before all the tests, read the config file and get a user token and
# save it to a file used by the main service script
class TestCoExpressionMethodsSetup(unittest.TestCase):
  def setUp(self):
    config = ConfigParser.RawConfigParser()
    config.read('ltest/test.cfg')
    user_id = config.get('CoExpressionTest','user')
    password = config.get('CoExpressionTest','password')
    token = Token(user_id=user_id, password=password)
    token_file = open('ltest/script_test/token.txt','w')
    token_file.write(token.token)

# Define all our other test cases here
class TestCoExpressionMethods(TestCoExpressionMethodsSetup): 

 def test_filter_genes(self):
        print("\n\n----------- test filter genes ----------")

        out =call(["run_CoExpression.sh",
        "ltest/script_test/test_view_heatmap_input.json",
        "ltest/script_test/test_view_heatmap_output.json",
        "ltest/script_test/token.txt"])

        # print error code of Implementation
        print(out);

        with open('ltest/script_test/test_view_heatmap_output.json') as o:
                output =json.load(o)
        pprint(output)

# def test_filter_genes(self):
#        print("\n\n----------- test filter genes ----------")
#
#        out =call(["run_CoExpression.sh",
#        "ltest/script_test/test_filter_genes_input.json",
#        "ltest/script_test/test_filter_genes_output.json",
#        "ltest/script_test/token.txt"])
#
#        # print error code of Implementation
#        print(out);
#
#        with open('ltest/script_test/test_filter_genes_output.json') as o:
#                output =json.load(o)
#        pprint(output)
#
# def test_coex_cluster(self):
#        print("\n\n----------- test constcoex_net_clust ----------")
#
#        out =call(["run_CoExpression.sh",
#        "ltest/script_test/test_coex_clust_input.json",
#        "ltest/script_test/test_coex_clust_output.json",
#        "ltest/script_test/token.txt"])
#
#        # print error code of Implementation
#        print(out);
#
#        with open('ltest/script_test/test_coex_clust_output.json') as o:
#                output =json.load(o)
#        pprint(output)


#start the tests if run as a script
if __name__ == '__main__':
    unittest.main()
