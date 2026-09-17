import unittest

import pandas as pd

from clustering_topology_layer import (
    annotate_topology_layer_evidence,
    build_topology_distance_matrix,
    extract_topology_signature,
    load_topology_layer_rules,
)


class TopologyLayerTests(unittest.TestCase):
    def test_rules_are_review_only(self):
        rules = load_topology_layer_rules()
        self.assertEqual(rules["mode"], "review_only")
        self.assertEqual(rules["rule_version"], "1.0.0-review-only")
        self.assertIn("Ring Profile", rules["topology_columns"])

    def test_signature_falls_back_to_ring_counts(self):
        row = pd.Series({
            "Ring Count": 2,
            "Aromatic Ring Count": 1,
        })
        signature = extract_topology_signature(row)
        self.assertEqual(signature["topology_class"], "Bicyclic")
        self.assertEqual(signature["ring_profile"], "2 ring(s), 1 aromatic ring(s)")

    def test_distance_matrix_prefers_shared_topology(self):
        df = pd.DataFrame([
            {
                "Topology Class": "Monocyclic",
                "Topology Modifiers": "Aromatic",
                "Ring System": "Monocyclic",
                "Ring Profile": "1 ring(s), 1 aromatic ring(s)",
                "Ring Count": 1,
                "Aromatic Ring Count": 1,
                "Alcohol Context": "Cyclic Alcohol",
            },
            {
                "Topology Class": "Monocyclic",
                "Topology Modifiers": "Aromatic",
                "Ring System": "Monocyclic",
                "Ring Profile": "1 ring(s), 1 aromatic ring(s)",
                "Ring Count": 1,
                "Aromatic Ring Count": 1,
                "Alcohol Context": "Cyclic Alcohol",
            },
            {
                "Topology Class": "Polycyclic",
                "Topology Modifiers": "Fused; Aromatic",
                "Ring System": "Polycyclic",
                "Ring Profile": "3 ring(s), 2 aromatic ring(s)",
                "Ring Count": 3,
                "Aromatic Ring Count": 2,
                "Alcohol Context": "Ring-substituted Alcohol",
            },
        ])
        matrix, signatures = build_topology_distance_matrix(df)
        self.assertEqual(len(signatures), 3)
        self.assertEqual(matrix[0, 1], 0.0)
        self.assertGreater(matrix[0, 2], 0.0)

    def test_cluster_annotation_reports_topology_cohesion(self):
        df = pd.DataFrame([
            {
                "Cluster ID": "Cluster_1",
                "Topology Class": "Monocyclic",
                "Topology Modifiers": "Aromatic",
                "Ring System": "Monocyclic",
                "Ring Profile": "1 ring(s), 1 aromatic ring(s)",
                "Ring Count": 1,
                "Aromatic Ring Count": 1,
                "Alcohol Context": "Cyclic Alcohol",
            },
            {
                "Cluster ID": "Cluster_1",
                "Topology Class": "Monocyclic",
                "Topology Modifiers": "Aromatic",
                "Ring System": "Monocyclic",
                "Ring Profile": "1 ring(s), 1 aromatic ring(s)",
                "Ring Count": 1,
                "Aromatic Ring Count": 1,
                "Alcohol Context": "Cyclic Alcohol",
            },
        ])
        result = annotate_topology_layer_evidence(df)
        self.assertEqual(result["Topology_Layer_Status"].iloc[0], "Calculated")
        self.assertIn(result["Topology_Layer_Cohesion"].iloc[0], {"Topologically coherent", "Topologically aligned"})
        self.assertIn("Mean pairwise topology distance", result["Topology_Layer_Reasons"].iloc[0])


if __name__ == "__main__":
    unittest.main()
