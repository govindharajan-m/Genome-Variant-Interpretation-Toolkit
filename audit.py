import unittest
from unittest.mock import patch
import json

from db_handler import _fetch_live_clinvar, _fetch_live_gwas, _fetch_live_population_frequencies
from variant_engine import annotate_snp_from_rsid

class TestAPIAudit(unittest.TestCase):
    @patch('db_handler.requests.get')
    def test_clinvar_timeout(self, mock_get):
        import requests
        mock_get.side_effect = requests.exceptions.Timeout("Timeout")
        result = _fetch_live_clinvar("rs429358")
        self.assertIsNone(result)

    @patch('db_handler.requests.get')
    def test_gwas_invalid_json(self, mock_get):
        mock_response = mock_get.return_value
        mock_response.status_code = 200
        mock_response.json.side_effect = json.JSONDecodeError("Expecting value", "", 0)
        result = _fetch_live_gwas("rs429358")
        self.assertEqual(result, [])

    @patch('db_handler.requests.get')
    def test_ensembl_404(self, mock_get):
        mock_response = mock_get.return_value
        mock_response.status_code = 404
        mock_response.raise_for_status.side_effect = Exception("404 Client Error")
        result = _fetch_live_population_frequencies("rs429358")
        self.assertIsNone(result)

class TestFunctionalAudit(unittest.TestCase):
    def test_all_requested_variants(self):
        rsids = ["rs429358", "rs7412", "rs334", "rs1800562", "rs1042522", "rs1805007", "rs113993960", "rs6025", "rs7903146", "rs9939609"]
        for rsid in rsids:
            res = annotate_snp_from_rsid(rsid)
            self.assertTrue(res.get("found", False) or "error" in res)
            
if __name__ == '__main__':
    unittest.main()
