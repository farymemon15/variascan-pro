import unittest
import os
from utils.profile_manager import load_profile, save_profile, PROFILE_FILE, DEFAULT_PROFILE


class TestProfileManager(unittest.TestCase):

    def setUp(self):
        self.orig_exists = PROFILE_FILE.exists()
        self.orig_content = PROFILE_FILE.read_text(encoding="utf-8") if self.orig_exists else None

    def tearDown(self):
        if self.orig_exists and self.orig_content:
            PROFILE_FILE.write_text(self.orig_content, encoding="utf-8")
        elif not self.orig_exists and PROFILE_FILE.exists():
            PROFILE_FILE.unlink()

    def test_load_default_profile(self):
        prof = load_profile()
        self.assertIn("name", prof)
        self.assertIn("email", prof)
        self.assertIn("upwork_active", prof)
        self.assertIn("linkedin_active", prof)
        self.assertIsInstance(prof["services"], list)

    def test_save_and_reload(self):
        prof = load_profile()
        prof["name"] = "Custom Bioinformatician"
        prof["linkedin_active"] = True
        prof["linkedin_url"] = "https://linkedin.com/in/test"
        
        saved = save_profile(prof)
        self.assertTrue(saved)
        
        reloaded = load_profile()
        self.assertEqual(reloaded["name"], "Custom Bioinformatician")
        self.assertTrue(reloaded["linkedin_active"])
        self.assertEqual(reloaded["linkedin_url"], "https://linkedin.com/in/test")


if __name__ == "__main__":
    unittest.main()
