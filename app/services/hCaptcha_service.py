import requests
from requests.exceptions import Timeout
import os

class HCaptchaService:
    hcaptcha_secret = os.environ.get("HCAPTCHA_SECRET")
    def __init__(self) -> None:
        pass

    def verify_captcha(self, captcha_token):
        try:
            response = requests.request('POST', f'''https://hcaptcha.com/siteverify''', timeout=2,
                                        data={
                                            'secret': self.hcaptcha_secret,
                                            'response': captcha_token,
                                        })
            try:
                assert response.status_code in [200, 204]
                print(result)
                result = response.json()
                if result['success']:
                    return True
                else:
                    return False
            except:
                return False
        except Timeout:
            return False