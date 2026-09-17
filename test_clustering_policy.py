import unittest
from pathlib import Path

from clustering_policy import load_clustering_policy


class ClusteringPolicyTests(unittest.TestCase):
    def test_policy_is_frozen_and_review_only(self):
        policy = load_clustering_policy()
        self.assertEqual(policy["mode"], "review_only")
        self.assertEqual(policy["rule_version"], "1.0.0-freeze")
        self.assertGreater(len(policy["distance_components"]), 0)
        total_weight = sum(component["weight"] for component in policy["distance_components"])
        self.assertAlmostEqual(total_weight, 1.0, places=6)

    def test_clustering_engine_source_references_the_policy_loader(self):
        source = Path("clustering_engine.py").read_text(encoding="utf-8")
        self.assertIn("load_clustering_policy", source)
        self.assertIn('df.attrs["clustering_policy"]', source)


if __name__ == "__main__":
    unittest.main()
