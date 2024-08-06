import logging
import os
import re
import requests
from requests import Timeout

from config import app_config

log = logging.getLogger(__name__)

SSL_VERIFY = True


def verify_captcha(captcha_token):
    try:
        response = requests.request('POST', f'''https://hcaptcha.com/siteverify''', timeout=2,
                                    data={
                                        'secret': app_config['auth']['hcaptcha_secret'],
                                        'response': captcha_token,
                                    }
                                    )
        try:
            assert response.status_code in [200, 204]
            result = response.json()
            if result['success']:
                return True
            else:
                log.error(f'''Invalid CAPTCHA''')
                return False
        except:
            log.error(f'''Could not verify CAPTCHA''')
            return False
    except Timeout:
        log.warning(f'''Could not verify CAPTCHA''')
        return False

def userinfo(access_token) -> dict:
    try:
        resp = requests.get(url=app_config['auth']['userInfoUrl'],
                            verify=SSL_VERIFY,
                            cookies={"_oauth2_proxy": access_token})
        resp.raise_for_status()
        user = resp.json()

        roles = []
        for grp in user['groups']:
            roles.append(grp)

        # sub = user['preferredUsername'].replace('@', '').replace('.', '')
        email = user['email']
        username = re.sub(r'[^a-zA-Z0-9]', '', email)

        # TODO: Hoping that oauth2-proxy enhances support for providing arbitrary token claims from OIDC
        # See https://github.com/oauth2-proxy/oauth2-proxy/issues/834
        return {
            'email': user['email'],
            'groups': roles,
            # 'family_name': user['family_name'],
            # 'given_name': user['given_name'],
            'sub': username
        }
    except Exception as e:
        log.debug(f'OAuth2 token verification failed for token={str(access_token)}')
        log.warning(f'OAuth2 token verification failed: {str(e)}')
        pass

def get_token_from_request_cookies(request):
    cookiename = app_config['auth']['cookiename']
    token = request.cookies.get(cookiename)
    #log.debug('Found get_token: ' + str(token))
    return token

# TODO: required scopes currently ignored
def validate_auth_cookie(request, required_scopes=[]):
    token = get_token_from_request_cookies(request)
    #log.debug('Found validate_token: ' + str(token))
    return userinfo(token)
